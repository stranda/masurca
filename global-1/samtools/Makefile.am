ACLOCAL_AMFLAGS = -I m4
bin_PROGRAMS = samtools
lib_LTLIBRARIES = libbam.la bcftools/libbcf.la

AM_CPPFLAGS = -I. -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE -D_CURSES_LIB=@CURSES_LIB@
AM_LDADD = -lm -lz -lpthread

library_includedir = $(includedir)/samtools
library_include_HEADERS = bam.h bam2bcf.h bam_endian.h bgzf.h		\
                          errmod.h faidx.h kaln.h khash.h klist.h	\
                          knetfile.h kprobaln.h kseq.h ksort.h		\
                          kstring.h razf.h sam.h sam_header.h		\
                          sample.h bam_tview.h

bcftools_includedir = $(library_includedir)/bcftools
bcftools_include_HEADERS = bcftools/bcf.h bcftools/kmin.h	\
                           bcftools/prob1.h

samtools_SOURCES = bam_tview.c bam_plcmd.c sam_view.c bam_rmdup.c	\
		   bam_rmdupse.c bam_mate.c bam_stat.c bam_color.c	\
		   bamtk.c kaln.c bam2bcf.c bam2bcf_indel.c errmod.c	\
		   sample.c cut_target.c phase.c bam2depth.c		\
		   padding.c bedcov.c bamshuf.c bam_tview_curses.c	\
		   bam_tview_html.c
samtools_LDADD = $(AM_LDADD) libbam.la bcftools/libbcf.la


libbam_la_SOURCES = bgzf.c kstring.c bam_aux.c bam.c bam_import.c	\
	            sam.c bam_index.c bam_pileup.c bam_lpileup.c	\
	            bam_md.c razf.c faidx.c bedidx.c knetfile.c		\
	            bam_sort.c sam_header.c bam_reheader.c kprobaln.c	\
	            bam_cat.c

bcftools_libbcf_la_SOURCES = bcftools/bcf.c bcftools/vcf.c		\
                             bcftools/bcfutils.c bcftools/prob1.c	\
                             bcftools/em.c bcftools/kfunc.c		\
                             bcftools/kmin.c bcftools/index.c		\
                             bcftools/fet.c bcftools/mut.c		\
                             bcftools/bcf2qcall.c bcftools/call1.c
