#!/usr/bin/env nextflow
/*
══════════════════════════════════════════════════════════════════════════════════
    STREAM — Streamlined Transcript Expression & RNA-seq Mapping
══════════════════════════════════════════════════════════════════════════════════
    Combined RNA-seq QC + Quantification Pipeline  •  Nextflow DSL2
    Ultra-minimalist design for solo bioinformaticians.

    Inspired by nf-core/rnaseq — fully self-contained in one file.

    Steps:
      0a. [Optional] GEO → SRR resolution  (NCBI E-utilities, Python stdlib)
      0b. [Optional] SRA download           (fasterq-dump)
      1.  FastQC — raw reads QC
      2.  fastp  — adapter trimming + quality filtering
      3.  FastQC — trimmed reads QC
      4.  [Optional] FastQ Screen — contamination screening
      5.  Sequence statistics               (seqtk)
      6.  [Optional] Kraken2 — taxonomic classification
      7.  [Optional] Salmon index build
      8.  [Optional] Salmon — transcript quantification
      9.  MultiQC — aggregate all reports

    GitHub: https://github.com/IPNP-BIPN/STREAM
══════════════════════════════════════════════════════════════════════════════════
*/

nextflow.enable.dsl = 2

// ════════════════════════════════════════════════════════════════════════════════
//  PROCESSES
// ════════════════════════════════════════════════════════════════════════════════

// ────────────────────────────────────────────────────────────────────────────────
//  PROCESS 0a : Résolution GEO (GSE/GSM) → accessions SRR
//  Utilise les NCBI E-utilities via Python stdlib (aucune dépendance externe).
//  Appelé uniquement lorsque --sra_ids contient des accessions GEO.
// ────────────────────────────────────────────────────────────────────────────────
process RESOLVE_GEO {
    tag "${geo_acc}"
    label 'process_low'
    maxForks 1

    input:
    val geo_acc

    output:
    path "srr_accessions.txt", emit: ids

    script:
    """
    #!/usr/bin/env python3
    import urllib.request, json, re, sys, time

    geo  = "${geo_acc}"
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    # Search SRA for the GEO accession
    url = f"{base}/esearch.fcgi?db=sra&term={geo}[Accession]&retmax=1000&retmode=json"
    with urllib.request.urlopen(url) as r:
        data = json.loads(r.read())
    ids = data.get("esearchresult", {}).get("idlist", [])
    if not ids:
        sys.exit(f"ERROR: No SRA entries found for {geo}")

    time.sleep(0.4)  # respect NCBI rate limit

    # Fetch run accessions via RunInfo CSV
    ids_str = ",".join(ids)
    url = f"{base}/efetch.fcgi?db=sra&id={ids_str}&rettype=runinfo&retmode=text"
    with urllib.request.urlopen(url) as r:
        text = r.read().decode()

    with open("srr_accessions.txt", "w") as fh:
        for line in text.strip().split("\\n")[1:]:
            cols = line.split(",")
            if cols and re.match(r"[SED]RR\\d+", cols[0]):
                fh.write(cols[0] + "\\n")
    """
}

// ────────────────────────────────────────────────────────────────────────────────
//  PROCESS 0b : Téléchargement FASTQ depuis SRA
//  Utilise prefetch + fasterq-dump (SRA Toolkit) pour récupérer les reads.
//  Gère automatiquement single-end et paired-end.
// ────────────────────────────────────────────────────────────────────────────────
process SRA_DOWNLOAD {
    tag "${accession}"
    label 'process_medium'
    publishDir "${params.outdir}/00_sra_fastq", mode: 'symlink'
    maxForks 3

    input:
    val accession

    output:
    tuple val(accession), path("*.fastq.gz"), emit: reads

    script:
    """
    # Prefetch .sra file (faster + more robust than streaming)
    prefetch ${accession} --max-size 100G

    # Convert to FASTQ (--split-files: PE → _1/_2, SE → single file)
    fasterq-dump ${accession} --threads ${task.cpus} --split-files --temp .

    # Compress
    pigz -p ${task.cpus} *.fastq 2>/dev/null || gzip *.fastq
    """
}

// ────────────────────────────────────────────────────────────────────────────────
//  PROCESS 1 : FastQC (brut)
//  Contrôle qualité des reads bruts avant tout traitement.
//  Les rapports HTML et ZIP sont collectés pour l'agrégation MultiQC.
// ────────────────────────────────────────────────────────────────────────────────
process FASTQC_RAW {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/01_fastqc_raw", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "*.html", emit: html
    path "*.zip",  emit: zip

    script:
    """
    fastqc --threads ${task.cpus} --quiet ${reads}
    """
}

// ────────────────────────────────────────────────────────────────────────────────
//  PROCESS 2 : fastp (trimming + QC)
//  Préprocesseur ultra-rapide all-in-one pour fichiers FASTQ :
//  contrôle qualité, trimming d'adaptateurs, filtrage qualité, pruning par read.
//  Pour les données paired-end, fastp détecte automatiquement les adaptateurs.
//  --qualified_quality_phred : score qualité minimum (Phred) pour qu'une base
//    soit considérée "qualifiée".
//  --length_required : élimine les reads plus courts que ce seuil après trimming.
//  --trim_poly_g / --trim_poly_x : supprime les queues poly-G (artefact NovaSeq)
//    et poly-X.
//  --cut_front / --cut_tail : trimming par fenêtre glissante en 5' et 3'.
//  Génère des rapports HTML et JSON détaillés avec métriques avant/après filtrage.
//  Ces rapports sont essentiels pour évaluer la qualité et l'efficacité du trimming.
// ────────────────────────────────────────────────────────────────────────────────
process FASTP {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/02_fastp", mode: 'copy', pattern: "*.{html,json}"
    publishDir "${params.outdir}/02_fastp/trimmed", mode: 'copy', pattern: "*_trimmed.fastq.gz", enabled: params.save_trimmed

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trimmed.fastq.gz"), emit: reads
    path "*.html",                                emit: html
    path "*.json",                                emit: json

    script:
    def prefix = meta.id
    if (meta.single_end) {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --out1 ${prefix}_R1_trimmed.fastq.gz \\
            --html ${prefix}_fastp.html \\
            --json ${prefix}_fastp.json \\
            --thread ${task.cpus} \\
            --qualified_quality_phred ${params.fastp_qualified_quality} \\
            --unqualified_percent_limit 10 \\
            --length_required ${params.fastp_length_required} \\
            --trim_poly_g \\
            --trim_poly_x \\
            --cut_front \\
            --cut_tail
        """
    } else {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${prefix}_R1_trimmed.fastq.gz \\
            --out2 ${prefix}_R2_trimmed.fastq.gz \\
            --html ${prefix}_fastp.html \\
            --json ${prefix}_fastp.json \\
            --thread ${task.cpus} \\
            --detect_adapter_for_pe \\
            --qualified_quality_phred ${params.fastp_qualified_quality} \\
            --unqualified_percent_limit 10 \\
            --length_required ${params.fastp_length_required} \\
            --trim_poly_g \\
            --trim_poly_x \\
            --cut_front \\
            --cut_tail
        """
    }
}

// ────────────────────────────────────────────────────────────────────────────────
//  PROCESS 3 : FastQC (nettoyé)
//  Contrôle qualité post-trimming pour vérifier l'amélioration de la qualité
//  et l'efficacité du trimming par fastp.
// ────────────────────────────────────────────────────────────────────────────────
process FASTQC_CLEAN {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/03_fastqc_clean", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "*.html", emit: html
    path "*.zip",  emit: zip

    script:
    """
    fastqc --threads ${task.cpus} --quiet ${reads}
    """
}

// ────────────────────────────────────────────────────────────────────────────────
//  PROCESS 4 : FastQ Screen (contamination)
//  Screening des reads trimmés contre plusieurs génomes de référence
//  (humain, souris, E. coli, adaptateurs…) pour détecter la contamination
//  croisée ou les contaminations exogènes.
//  Utilise bowtie2 comme aligneur. Le paramètre subset_size permet un test rapide.
// ────────────────────────────────────────────────────────────────────────────────
process FASTQ_SCREEN {
    tag "${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/04_fastq_screen", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "*_screen.*", emit: report

    script:
    def subset_opt = params.subset_size > 0 ? "--subset ${params.subset_size}" : ""
    """
    fastq_screen \\
        --conf ${params.fastq_screen_conf} \\
        --aligner bowtie2 \\
        --threads ${task.cpus} \\
        ${subset_opt} \\
        ${reads[0]}
    """
}

// ────────────────────────────────────────────────────────────────────────────────
//  PROCESS 5 : Statistiques de séquences (seqtk)
//  Génération de statistiques globales pour chaque échantillon :
//  nombre de reads, bases totales, longueur moyenne, composition en bases.
//  Utilise seqtk fqchk pour un profil rapide et complet des fichiers FASTQ.
// ────────────────────────────────────────────────────────────────────────────────
process SEQTK_STATS {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/05_statistics", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "*.fqchk.txt", emit: stats

    script:
    """
    for fq in ${reads}; do
        base=\$(basename \$fq .fastq.gz)
        seqtk fqchk \$fq > \${base}.fqchk.txt
    done
    """
}

// ────────────────────────────────────────────────────────────────────────────────
//  PROCESS 6 : Kraken2 — classification taxonomique (optionnel)
//  Identification des espèces présentes dans les reads trimmés via une approche
//  k-mer sur une base de données Kraken2 (standard ou custom).
//  Les rapports (.report) sont compatibles MultiQC.
//  Nécessite une base Kraken2 complète (taxo.k2d, hash.k2d, opts.k2d).
// ────────────────────────────────────────────────────────────────────────────────
process KRAKEN2 {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/06_kraken2", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path "*.report", emit: report
    path "*.out",    emit: output

    script:
    def prefix = meta.id
    if (meta.single_end) {
        """
        kraken2 \\
            --db ${params.kraken2_db} \\
            --threads ${task.cpus} \\
            ${reads[0]} \\
            --report ${prefix}.report \\
            --output ${prefix}.out
        """
    } else {
        """
        kraken2 \\
            --db ${params.kraken2_db} \\
            --threads ${task.cpus} \\
            --paired ${reads[0]} ${reads[1]} \\
            --report ${prefix}.report \\
            --output ${prefix}.out
        """
    }
}

// ────────────────────────────────────────────────────────────────────────────────
//  PROCESS 7a : Téléchargement du transcriptome Ensembl
//  Ensembl fournit une annotation complète des gènes et transcrits humains.
//  Le fichier cDNA FASTA contient toutes les séquences de transcrits connues
//  pour l'assemblage spécifié (par défaut GRCh38).
//  Ce fichier est nécessaire pour construire l'index Salmon.
//  Le script vérifie si le fichier existe déjà (storeDir) pour éviter
//  les téléchargements redondants.
//  NOTE : URL hardcodée pour Homo sapiens. Pour d'autres organismes,
//  fournir --transcriptome_fasta directement.
// ────────────────────────────────────────────────────────────────────────────────
process DOWNLOAD_TRANSCRIPTOME {
    label 'process_low'
    storeDir "${params.outdir}/reference"

    output:
    path "transcriptome.fa", emit: fasta

    script:
    def assembly = params.genome
    def release  = params.ensembl_release
    """
    curl -L "https://ftp.ensembl.org/pub/release-${release}/fasta/homo_sapiens/cdna/Homo_sapiens.${assembly}.cdna.all.fa.gz" \\
        -o transcriptome.fa.gz
    gunzip transcriptome.fa.gz
    """
}

// ────────────────────────────────────────────────────────────────────────────────
//  PROCESS 7b : Construction de l'index Salmon
//  Salmon utilise un index k-mer spécialisé pour le quasi-mapping des reads
//  vers les transcrits. Cet index accélère considérablement la quantification
//  par rapport à l'alignement traditionnel.
//  L'index est construit une seule fois et réutilisé pour tous les échantillons.
//  Utilise storeDir pour persister l'index entre les exécutions du pipeline.
// ────────────────────────────────────────────────────────────────────────────────
process SALMON_INDEX {
    label 'process_high'
    storeDir "${params.outdir}/reference"

    input:
    path fasta

    output:
    path "salmon_index", emit: index

    script:
    """
    salmon index \\
        -t ${fasta} \\
        -i salmon_index \\
        -p ${task.cpus}
    """
}

// ────────────────────────────────────────────────────────────────────────────────
//  PROCESS 8 : Quantification Salmon
//  Salmon effectue le quasi-mapping et la quantification de l'abondance
//  des transcrits à partir des reads RNA-seq.
//  Contrairement aux outils d'alignement traditionnels, Salmon utilise une
//  approche de mapping léger, à la fois plus rapide et plus précise.
//  -l A : détection automatique du type de librairie (ISR, ISF, IU…).
//  --gcBias : correction des biais GC au niveau des fragments.
//  --validateMappings : validation supplémentaire des mappings via scoring
//    d'alignement, améliore la précision au coût d'un runtime légèrement plus long.
//  Résultats : abondances en TPM et comptages estimés (quant.sf),
//  directement utilisables pour l'analyse d'expression différentielle (DESeq2, etc.).
// ────────────────────────────────────────────────────────────────────────────────
process SALMON_QUANT {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/07_salmon", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}"), emit: results

    script:
    def prefix = meta.id
    if (meta.single_end) {
        """
        salmon quant \\
            -i ${index} \\
            -l A \\
            -r ${reads[0]} \\
            -p ${task.cpus} \\
            -o ${prefix} \\
            --gcBias \\
            --validateMappings
        """
    } else {
        """
        salmon quant \\
            -i ${index} \\
            -l A \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -p ${task.cpus} \\
            -o ${prefix} \\
            --gcBias \\
            --validateMappings
        """
    }
}

// ────────────────────────────────────────────────────────────────────────────────
//  PROCESS 9 : MultiQC — agrégation de tous les rapports
//  MultiQC agrège les rapports de FastQC, fastp, Salmon, Kraken2 et FastQ Screen
//  en un seul rapport HTML interactif.
//  Rapport principal : ${outdir}/08_multiqc/multiqc_report.html
// ────────────────────────────────────────────────────────────────────────────────
process MULTIQC {
    label 'process_low'
    publishDir "${params.outdir}/08_multiqc", mode: 'copy'

    input:
    path('*')

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data

    script:
    """
    multiqc . --force
    """
}


// ════════════════════════════════════════════════════════════════════════════════
//  WORKFLOW
// ════════════════════════════════════════════════════════════════════════════════

workflow {

    // ── Validate inputs ─────────────────────────────────────────────────────
    if (!params.input && !params.fastq_dir && !params.sra_ids) {
        error """
        ╔══════════════════════════════════════════════════════════════════╗
        ║  STREAM requires one of:                                        ║
        ║    --input      samplesheet.csv                                 ║
        ║    --fastq_dir  /path/to/fastqs                                 ║
        ║    --sra_ids    SRR1234,SRR5678  (or file with accessions)      ║
        ╚══════════════════════════════════════════════════════════════════╝
        """.stripIndent()
    }

    // ── Banner ──────────────────────────────────────────────────────────────
    log.info """
    ╔══════════════════════════════════════════════════════════════════╗
    ║   S T R E A M                                                    ║
    ║   Streamlined Transcript Expression & RNA-seq Mapping            ║
    ║   v${workflow.manifest.version ?: '1.0.0'}                                                       ║
    ╚══════════════════════════════════════════════════════════════════╝
    Input        : ${params.input ?: params.fastq_dir ?: params.sra_ids}
    Output       : ${params.outdir}
    Salmon       : ${params.run_salmon}
    FastQ Screen : ${params.run_fastq_screen}
    Kraken2      : ${params.run_kraken2}
    CPUs         : ${params.max_cpus}
    ──────────────────────────────────────────────────────────────────
    """.stripIndent()

    // ════════════════════════════════════════════════════════════════════════
    //  INPUT CHANNEL ASSEMBLY
    // ════════════════════════════════════════════════════════════════════════

    ch_raw_reads = Channel.empty()

    // ── Mode SRA/GEO ────────────────────────────────────────────────────────
    if (params.sra_ids) {

        // Parse accessions: comma-separated string OR file path (one per line)
        if (file(params.sra_ids).exists()) {
            ch_accessions = Channel.fromPath(params.sra_ids)
                .splitText()
                .map { it.trim() }
                .filter { it && !it.startsWith('#') }
        } else {
            ch_accessions = Channel
                .from(params.sra_ids.toString().split(',').collect { it.trim() })
        }

        // Separate GEO accessions (GSE/GSM) from direct SRA accessions (SRR/ERR/DRR)
        ch_geo        = ch_accessions.filter { it ==~ /^GS[EM]\d+$/ }
        ch_sra_direct = ch_accessions.filter { it ==~ /^[SED]RR\d+$/ }

        // Resolve GEO → SRR accessions via NCBI E-utilities
        RESOLVE_GEO(ch_geo)
        ch_srr_from_geo = RESOLVE_GEO.out.ids
            .splitText()
            .map { it.trim() }
            .filter { it }

        // Merge all SRR accessions and download
        ch_all_sra = ch_sra_direct.mix(ch_srr_from_geo)
        SRA_DOWNLOAD(ch_all_sra)

        // Build [meta, reads] tuples — auto-detect PE/SE from file count
        ch_raw_reads = SRA_DOWNLOAD.out.reads
            .map { accession, files ->
                def flist  = files instanceof List ? files : [files]
                def sorted = flist.sort { it.name }
                // Filter orphan reads (_3.fastq.gz) if present
                sorted = sorted.findAll { !(it.name =~ /_3\.fastq/) }
                def meta = [id: accession, single_end: sorted.size() == 1]
                [meta, sorted]
            }

    // ── Mode fastq_dir (auto-détection PE/SE) ──────────────────────────────
    } else if (params.fastq_dir) {

        ch_raw_reads = Channel
            .fromFilePairs("${params.fastq_dir}/*_R{1,2}_001.fastq.gz", size: -1)
            .map { sample, files ->
                def meta = [id: sample, single_end: files.size() == 1]
                [meta, files.sort { it.name }]
            }

    // ── Mode samplesheet CSV ────────────────────────────────────────────────
    //  Format: sample,fastq_1,fastq_2   (fastq_2 vide = single-end)
    } else if (params.input) {

        ch_raw_reads = Channel
            .fromPath(params.input)
            .splitCsv(header: true, strip: true)
            .map { row ->
                def has_r2 = row.fastq_2?.trim()
                def meta   = [id: row.sample, single_end: !has_r2]
                def reads  = has_r2
                    ? [file(row.fastq_1), file(row.fastq_2)]
                    : [file(row.fastq_1)]
                [meta, reads]
            }
    }

    // ════════════════════════════════════════════════════════════════════════
    //  PIPELINE EXECUTION
    // ════════════════════════════════════════════════════════════════════════

    // Collect all QC outputs for MultiQC aggregation
    ch_multiqc = Channel.empty()

    // ── 1. FastQC raw ───────────────────────────────────────────────────────
    if (!params.skip_fastqc) {
        FASTQC_RAW(ch_raw_reads)
        ch_multiqc = ch_multiqc.mix(FASTQC_RAW.out.zip)
    }

    // ── 2. fastp (always runs — core trimming step) ─────────────────────────
    FASTP(ch_raw_reads)
    ch_trimmed = FASTP.out.reads
    ch_multiqc = ch_multiqc.mix(FASTP.out.json)

    // ── 3. FastQC clean ─────────────────────────────────────────────────────
    if (!params.skip_fastqc) {
        FASTQC_CLEAN(ch_trimmed)
        ch_multiqc = ch_multiqc.mix(FASTQC_CLEAN.out.zip)
    }

    // ── 4. FastQ Screen (optional — needs config) ───────────────────────────
    if (params.run_fastq_screen && params.fastq_screen_conf) {
        FASTQ_SCREEN(ch_trimmed)
        ch_multiqc = ch_multiqc.mix(FASTQ_SCREEN.out.report)
    }

    // ── 5. Sequence statistics ──────────────────────────────────────────────
    SEQTK_STATS(ch_trimmed)

    // ── 6. Kraken2 (optional — needs database) ─────────────────────────────
    if (params.run_kraken2 && params.kraken2_db) {
        KRAKEN2(ch_trimmed)
        ch_multiqc = ch_multiqc.mix(KRAKEN2.out.report)
    }

    // ── 7–8. Salmon quantification (optional) ──────────────────────────────
    if (params.run_salmon) {

        // Prepare Salmon index (3 paths: pre-built > local fasta > download)
        if (params.salmon_index) {
            ch_index = Channel.value(file(params.salmon_index))
        } else {
            if (params.transcriptome_fasta) {
                ch_fasta = Channel.value(file(params.transcriptome_fasta))
            } else {
                DOWNLOAD_TRANSCRIPTOME()
                ch_fasta = DOWNLOAD_TRANSCRIPTOME.out.fasta
            }
            SALMON_INDEX(ch_fasta)
            ch_index = SALMON_INDEX.out.index
        }

        SALMON_QUANT(ch_trimmed, ch_index.first())
        ch_multiqc = ch_multiqc.mix(
            SALMON_QUANT.out.results.map { meta, dir -> dir }
        )
    }

    // ── 9. MultiQC aggregation ──────────────────────────────────────────────
    MULTIQC(ch_multiqc.collect())
}


// ════════════════════════════════════════════════════════════════════════════════
//  ON COMPLETE
// ════════════════════════════════════════════════════════════════════════════════

workflow.onComplete {
    log.info """
    ══════════════════════════════════════════════════════════════════
      STREAM complete !
      Status   : ${workflow.success ? '✅ SUCCESS' : '❌ FAILED'}
      Duration : ${workflow.duration}
      Output   : ${params.outdir}
      Report   : ${params.outdir}/08_multiqc/multiqc_report.html
    ══════════════════════════════════════════════════════════════════
    """.stripIndent()
}


/*
 * ════════════════════════════════════════════════════════════════════════════════
 *  EXEMPLES D'EXÉCUTION
 * ════════════════════════════════════════════════════════════════════════════════
 *
 *  # Depuis un dossier de FASTQs (auto-détection PE/SE) :
 *  nextflow run main.nf --fastq_dir /path/to/fastqs --outdir results -resume
 *
 *  # Depuis un samplesheet CSV :
 *  nextflow run main.nf --input samplesheet.csv --outdir results -resume
 *
 *  # Depuis des accessions SRA :
 *  nextflow run main.nf --sra_ids "SRR1234567,SRR1234568" --outdir results -resume
 *
 *  # Depuis un numéro GEO (GSE) :
 *  nextflow run main.nf --sra_ids GSE123456 --outdir results -resume
 *
 *  # Avec un index Salmon pré-construit :
 *  nextflow run main.nf \
 *      --fastq_dir fastqs \
 *      --salmon_index /path/to/salmon_index \
 *      --outdir results -resume
 *
 *  # Pipeline complet avec toutes les options QC :
 *  nextflow run main.nf \
 *      --fastq_dir /path/to/fastqs \
 *      --run_salmon \
 *      --run_fastq_screen --fastq_screen_conf /path/to/fastq_screen.conf \
 *      --run_kraken2 --kraken2_db /path/to/kraken2_db \
 *      --outdir results \
 *      -resume
 *
 *  # Depuis GitHub :
 *  nextflow run IPNP-BIPN/STREAM --fastq_dir /path/to/fastqs --outdir results -resume
 *
 * ════════════════════════════════════════════════════════════════════════════════
 */
