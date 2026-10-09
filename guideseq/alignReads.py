"""
alignReads
"""

import subprocess
import os
import logging

logger = logging.getLogger('root')
logger.propagate = False


def alignReads(HG19_path, read1, read2, outfile,njobs=6,umi_tools="umi_tools",samtools="samtools",bwa="bwa"):

	sample_name = os.path.basename(outfile).split('.')[0]
	output_folder = os.path.dirname(outfile)
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)

	sample_alignment_paths = {}

	# When several samples share one untreated control library there is no reason to
	# align it more than once, and an interrupted run should be resumable. Skip the
	# work if the deduplicated alignment for this name is already present.
	final_sam = os.path.join(output_folder, sample_name + '.dedup.sam')
	if os.path.isfile(final_sam) and os.path.getsize(final_sam) > 0:
		logger.info('Alignment already present for %s, skipping (%s)', sample_name, final_sam)
		return

	# Check if genome is already indexed by bwa
	index_files_extensions = ['.pac', '.amb', '.ann', '.bwt', '.sa']

	genome_indexed = True
	for extension in index_files_extensions:
		if not os.path.isfile(HG19_path + extension):
			genome_indexed = False
			break

	# If the genome is not already indexed, index it
	if not genome_indexed:
		logger.info('Genome index files not detected. Running BWA to generate indices.')
		bwa_index_command = '{0} index {1}'.format(bwa, HG19_path)
		logger.info('Running bwa command: %s', bwa_index_command)
		subprocess.call(bwa_index_command.split())
		logger.info('BWA genome index generated')
	else:
		logger.info('BWA genome index found.')

	# Run paired end alignment against the genome
	logger.info('Running paired end mapping for {0}'.format(sample_name))
	# Stream bwa output straight into coordinate sorting. Avoids writing the
	# intermediate .sam/.bam, which for a full-depth GUIDE-seq library is tens of
	# gigabytes; the sorted BAM produced is identical to the previous three-step form.
	bwa_alignment_command = '{0} mem -t {1} {2} {3} {4} | {5} sort -@ 2 -m 1G -o {6}/{7}.st.bam -'.format(bwa,
														 njobs,
														 HG19_path,
														 read1,
														 read2,
														 samtools,
														 output_folder,
														 sample_name
														 )

	logger.info(bwa_alignment_command)
	if os.system(bwa_alignment_command) != 0:
		raise RuntimeError('bwa/samtools sort failed for %s' % sample_name)
	command = "{0} index {1}/{2}.st.bam".format(samtools,output_folder,sample_name)
	logger.info(command)
	os.system(command)
	# command = "{0} dedup --stdin={1}/{2}.st.bam --log={1}/{2}.dedup.log --output-stats={1}/{2}.stats.tsv --paired > {1}/{2}.dedup.bam;{3} index {1}/{2}.dedup.bam;{3} view {1}/{2}.dedup.bam > {1}/{2}.dedup.sam".format(umi_tools,output_folder,sample_name,samtools)
	# removing --output-stats, for large memory usage
	command = "{0} dedup --method unique --stdin={1}/{2}.st.bam --log={1}/{2}.dedup.log --paired > {1}/{2}.dedup.bam;{3} index {1}/{2}.dedup.bam;{3} view {1}/{2}.dedup.bam > {1}/{2}.dedup.sam".format(umi_tools,output_folder,sample_name,samtools)
	logger.info(command)
	os.system(command)
	logger.info('Paired end mapping for {0} completed.'.format(sample_name))
