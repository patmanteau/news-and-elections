import shutil

assert len(snakemake.input) == 1
assert len(snakemake.output) == 1

shutil.copy(snakemake.input[0], snakemake.output[0])
