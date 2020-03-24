## run under /av/NGSbackups
# NAPdemux --table demultiplexing_tables/LB3806975_RNMB_demux.tsv --output ../QMRmeta/mbc/00_demux_LB/ --table demultiplexing_tables/LB3806975_RNBC_demux.tsv --convert demultiplexing_tables/LB3806975_convert.csv --pairs strict --verbose AVlab_2018-05_HiSeq_LB3806975_amplicons/*.fastq

# NAPdemux --table demultiplexing_tables/NHMMar2019_RNMB_demux.tsv --output ../QMRmeta/mbc/00_demux_NHM/ --table demultiplexing_tables/NHMMar2019_RNBC_demux.tsv --convert demultiplexing_tables/NHMMar2019_convert.tsv --pairs strict --verbose AVlab_2019-03_MiSeq_NHMMar2019/*.fastq
NAPdemux --table demultiplexing_tables/LB3806975_RNMB_demux.tsv --output ../QMRmeta/mbc/00_demux_LB/ --table demultiplexing_tables/LB3806975_RNBC_demux.tsv --convert demultiplexing_tables/LB3806975_convert.csv --pairs strict --verbose LB3806975_amplicons/*.fastq |& tee -a ../QMRmeta/output/00_demux_LB.log

NAPdemux --table demultiplexing_tables/NHMMar2019_RNMB_demux.tsv --output ../QMRmeta/mbc/00_demux_NHM/ --table demultiplexing_tables/NHMMar2019_RNBC_demux.tsv --convert demultiplexing_tables/NHMMar2019_convert.tsv --pairs strict --verbose NHMMar2019/*.fastq |& tee -a ../QMRmeta/output/00_demux_NHM.log





# without remove and write singletons.
NAPdemux --table demultiplexing_tables/LB3806975_RNMB_demux.tsv --output ../QMRmeta/mbc/00_demux_LB/ --table demultiplexing_tables/LB3806975_RNBC_demux.tsv --convert demultiplexing_tables/LB3806975_convert.csv --verbose LB3806975_amplicons/*.fastq |& tee -a ../QMRmeta/output/00_demux_LB.log

NAPdemux --table demultiplexing_tables/NHMMar2019_RNMB_demux.tsv --output ../QMRmeta/mbc/00_demux_NHM/ --table demultiplexing_tables/NHMMar2019_RNBC_demux.tsv --convert demultiplexing_tables/NHMMar2019_convert.tsv --verbose NHMMar2019/*.fastq |& tee -a ../QMRmeta/output/00_demux_NHM.log


# with error tolerance = 2
NAPdemux --table demultiplexing_tables/LB3806975_RNMB_demux.tsv --output ../QMRmeta/mbc/00_demux_LB/ --table demultiplexing_tables/LB3806975_RNBC_demux.tsv --convert demultiplexing_tables/LB3806975_convert.csv --verbose --errortolerance 2 LB3806975_amplicons/*.fastq |& tee -a ../QMRmeta/output/00_demux_LB.log

NAPdemux --table demultiplexing_tables/NHMMar2019_RNMB_demux.tsv --output ../QMRmeta/mbc/00_demux_NHM/ --table demultiplexing_tables/NHMMar2019_RNBC_demux.tsv --convert demultiplexing_tables/NHMMar2019_convert.tsv --errortolerance 2 --verbose NHMMar2019/*.fastq |& tee -a ../QMRmeta/output/00_demux_NHM.log


## demux missing data:
NAPdemux --table demultiplexing_tables/LB3806975_RNMB_demux.tsv --output ../QMRmeta/mbc/00_demux_LB/ --table demultiplexing_tables/LB3806975_RNBC_demux.tsv --convert demultiplexing_tables/LB3806975_convert.csv --pairs strict --verbose LB3806975_amplicons/*.fastq |& tee -a ../QMRmeta/output/00_demux_LB.log

NAPdemux --table demultiplexing_tables/NHMMar2019_RNMB_demux.tsv --output ../QMRmeta/mbc/00_demux_NHM/ --table demultiplexing_tables/NHMMar2019_RNBC_demux.tsv --convert demultiplexing_tables/NHMMar2019_convert.tsv --pairs strict --verbose NHMMar2019/*.fastq |& tee -a ../QMRmeta/output/00_demux_NHM.log