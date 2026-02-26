from viral_usher import viral_usher_build
from viral_usher import ncbi_helper
from viral_usher import config
ncbi = ncbi_helper.NcbiHelper()
config_contents = config.parse_config("config.toml")
ref_acc, ref_fasta, ref_gbff, ref_length, ref_segment = viral_usher_build.get_reference(config_contents, ncbi)
