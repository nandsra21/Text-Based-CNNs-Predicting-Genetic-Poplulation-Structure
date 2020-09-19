"""
Rules to build simulated datasets.
"""

#
# Rules for simulated datasets
#

EMBEDDING_METHODS = [
    "pca",
    "mds",
    "t-sne",
    "umap"
]


rule all:
    input:
        expand("../../../../d/santa-sim/simulated_data/results/scatterplot_{method}.png", method=EMBEDDING_METHODS),
        expand("../../../../d/santa-sim/simulated_data/results/KDEDensity_{method}.png", method=EMBEDDING_METHODS),
        "../../../../d/santa-sim/simulated_data/results/tree_raw.nwk"

rule run_simulation:
    input:
        simulation_config = "least_recombinant_h3n2ha.xml"
    output:
        sequences = "../../../../d/santa-sim/simulated_data/simulated_HA_sequences.fasta"
    
    shell:
        """
        java -jar ../dist/santa.jar {input.simulation_config}
        """


rule parse_simulated_sequences:
    input:
        sequences = rules.run_simulation.output.sequences
    output:
        sequences = "../../../../d/santa-sim/simulated_data/sequences.fasta",
        metadata = "../../../../d/santa-sim/simulated_data/metadata.tsv"
    params:
        fasta_fields = "strain generation fitness"
    
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = "../../../../d/santa-sim/simulated_data/sequences.fasta"
    output:
        tree = "../../../../d/santa-sim/simulated_data/results/tree_raw.nwk"
    
    threads: 4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads {threads}
        """

rule create_distance_matrix:
    message: "creating the distance matrix to be used in the rest of the analysis"
    input:
        alignment = "../../../../d/santa-sim/simulated_data/sequences.fasta"
    output:
        output = "../../../../d/santa-sim/simulated_data/results/distance_matrix.csv"
    
    shell:
        """
        python3 scripts/hamming_distance_from_fasta.py \
            --alignment {input.alignment} \
            --output {output.output}
        """

rule embed_pca:
    message: "Creating the embedding (dataframe, node JSON) for PCA"
    input:
        alignment = "../../../../d/santa-sim/simulated_data/sequences.fasta"
    output:
        dataframe = "../../../../d/santa-sim/simulated_data/results/embed_pca.csv"
    params:
        components = 10
    
    shell:
        """
        python3 scripts/embed.py \
            --alignment {input.alignment} \
            --cluster True \
            --output-node-data {output.node_data} \
            --output-dataframe {output.dataframe} \
            pca \
            --components {params.components}
        """

rule embed_tsne:
    message: "Creating the embedding (dataframe, node JSON) for t-SNE"
    input:
        distance_matrix = rules.create_distance_matrix.output.output
    output:
        dataframe = "../../../../d/santa-sim/simulated_data/results/embed_t-sne.csv"
    params:
        perplexity = 15,
        learning_rate = 100
    
    shell:
        """
        python3 scripts/embed.py \
            --distance-matrix {input.distance_matrix} \
            --cluster True \
            --output-dataframe {output.dataframe} \
            t-sne \
            --perplexity {params.perplexity} \
            --learning-rate {params.learning_rate}
        """

rule embed_umap:
    message: "Creating the embedding (dataframe, node JSON) for UMAP"
    input:
        distance_matrix = rules.create_distance_matrix.output.output
    output:
        dataframe = "../../../../d/santa-sim/simulated_data/results/embed_umap.csv"
    params:
        nearest_neighbors = 200,
        min_dist = .05
    
    shell:
        """
        python3 scripts/embed.py \
            --distance-matrix {input.distance_matrix} \
            --cluster True \
            --output-dataframe {output.dataframe} \
            umap \
            --nearest-neighbors {params.nearest_neighbors} \
            --min-dist {params.min_dist}
        """

rule embed_mds:
    message: "Creating the embedding (dataframe, node JSON) for MDS"
    input:
        distance_matrix = rules.create_distance_matrix.output.output
    output:
        dataframe = "../../../../d/santa-sim/simulated_data/results/embed_mds.csv"
    params:
        components = 10
    
    shell:
        """
        python3 scripts/embed.py \
            --distance-matrix {input.distance_matrix} \
            --cluster True \
            --output-dataframe {output.dataframe} \
            mds \
            --components {params.components} \
        """

rule scatterplot:
    message: "Creating the scatterplot (PNG, dataframe)"
    input:
        distance_matrix = "../../../../d/santa-sim/simulated_data/results/distance_matrix.csv",
        embedding = "../../../../d/santa-sim/simulated_data/results/embed_{method}.csv"
    output:
        figure = "../../../../d/santa-sim/simulated_data/results/scatterplot_{method}.png",
    
    shell:
        """
        python3 scripts/scatterplot.py \
            --distance {input.distance_matrix} \
            --embedding {input.embedding} \
            --method {wildcards.method} \
            --output-figure {output.figure}
        """

def _get_embedding_columns_by_wildcards(wildcards):
    method = wildcards.method.replace("-", "")

    if method in ("pca", "mds"):
        return f"{method}1 {method}2"
    else:
        return f"{method}_x {method}_y"

rule KDE_density:
    message: "creating the KDE density plot"
    input:
        embedding = "../../../../d/santa-sim/simulated_data/results/embed_{method}.csv",
        metadata = "../../../../d/santa-sim/simulated_data/metadata.tsv"
    output:
        figure = "../../../../d/santa-sim/simulated_data/results/KDEDensity_{method}.png",
        dataframe = "../../../../d/santa-sim/simulated_data/results/KDEDensity_{method}.csv"
    params:
        embedding_columns = _get_embedding_columns_by_wildcards,
        differentiator_column = "labels"
    
    shell:
        """
        python3 scripts/within_vs_between_status.py \
            --embedding {input.embedding} \
            --metadata {input.metadata} \
            --method {wildcards.method} \
            --embedding-columns {params.embedding_columns} \
            --differentiator-column {params.differentiator_column} \
            --output-figure {output.figure} \
            --output-dataframe {output.dataframe}
        """
rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    
    shell:
        "rm -rfv {params}"
