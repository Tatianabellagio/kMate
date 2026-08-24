#!/usr/bin/env Rscript
# Combined CAM5 3'-REGION figure: LFMM GEA Manhattan (top, frameless) over the
# 17-haplotype 3' pangenome panel (bottom), sharing the genomic x-axis. Dashed
# droplines connect the Bonferroni-significant LFMM variants down into the
# haplotype map so you can see which haplotypes carry the significant sites.
suppressMessages({library(gggenomes); library(dplyr); library(readr)
                  library(patchwork); library(scales)})
HERE <- "/global/scratch/users/tbellg/kmate/analysis/grenenet_selection/r2_gea_nonsnp/cam5_replication/pang_cam5"
infct <- function(x) factor(x, levels = unique(x))
WIN_LO <- 11531800
HAP_LO <- 11533904-WIN_LO; HAP_HI <- 11534143-WIN_LO
SQUARES <- length(commandArgs(TRUE))>0 && commandArgs(TRUE)[1]=="squares"
SUF <- if (SQUARES) "_squares" else ""

reg <- read_tsv(file.path(HERE,"gg_region.tsv"), show_col_types=FALSE)
RLO <- reg$x[reg$name=="region_lo"]; RHI <- reg$x[reg$name=="region_hi"]
B1MEAN <- reg$x[reg$name=="bio1_mean"]

# significant LFMM variants (Bonferroni over the gene's 71 tested variants)
av <- read_tsv(file.path(HERE,"cam5_variants_gen9.tsv"), show_col_types=FALSE) |> filter(!is.na(lfmm_p))
bonf <- 0.05/nrow(av)
sigx <- (av |> filter(lfmm_p < bonf) |> mutate(x=pos-WIN_LO))$x      # significant gene-x positions
dl <- function() geom_vline(xintercept=sigx, linetype="dashed", colour="grey45",
                            linewidth=0.45, alpha=0.55)

## ---------------- bottom: 3' region haplotype panel ----------------
seqs <- read_tsv(file.path(HERE,"gg_rhapseqs.tsv"),  show_col_types=FALSE) |>
  mutate(seq_id=as.character(seq_id)) |> filter(kind=="hap") |> mutate(bin_id=infct(seq_id))
model <- read_tsv(file.path(HERE,"gg_model.tsv"),    show_col_types=FALSE) |> mutate(seq_id=as.character(seq_id))
vars  <- read_tsv(file.path(HERE,"gg_rhapfeats.tsv"), show_col_types=FALSE) |> mutate(seq_id=as.character(seq_id))

p  <- gggenomes(seqs=seqs); ly <- pull_seqs(p); hap <- ly
ymin <- min(ly$y); ymax <- max(ly$y)
v  <- vars |> left_join(select(ly, seq_id, y), by="seq_id") |> filter(!is.na(y))
vp <- v |> filter(vtype %in% c("snp","mnp")) |> mutate(xm=(start+end)/2)
vi <- v |> filter(vtype=="ins") |> mutate(yo=y+0.30)
vd <- v |> filter(vtype=="del") |> mutate(xe=pmax(end, start+1))
dlab <- vd |> filter(size>=5) |> distinct(start, xe, size) |> mutate(xm=(start+xe)/2)
ISO <- "AT2G27030.3"; gyv <- ymax + 1.7
ext <- model |> filter(seq_id==ISO) |> summarise(x0=min(start), x1=max(end)) |> mutate(gy=gyv)
mm  <- model |> filter(seq_id==ISO) |> mutate(gy=gyv)
glab_top <- ymax + 2.9
BARGAP <- 10; BARMAX <- 75; K <- BARMAX/sqrt(max(hap$n_eco))
hap <- hap |> mutate(xb1=RLO-BARGAP, xb0=RLO-BARGAP - sqrt(n_eco)*K)
barleft <- RLO - BARGAP - BARMAX - 20; xleft <- barleft
ref_id <- setdiff(seqs$seq_id, unique(vars$seq_id))[1]
ref_y  <- ly$y[ly$seq_id==ref_id]; ref_n <- hap$n_eco[hap$seq_id==ref_id]
ref_xc <- with(filter(hap, seq_id==ref_id), (xb0+xb1)/2)
brk <- c(2000,2100,2200,2300,2400,2500)

# variant glyph layers + colour scale depend on mode
if (SQUARES) {
  vv <- v |> mutate(grp=ifelse(vtype=="snp","SNP","indel"), xm=(start+end)/2)
  var_layers <- list(geom_point(aes(x=xm, y=y, colour=grp, shape=grp), data=vv, size=2.2),
                     scale_shape_manual(values=c(SNP=16, indel=15), name="variant"))
  col_scale <- scale_colour_manual(values=c(SNP="grey30", indel="grey55", CDS="#2b2b2b", UTR="#9e9e9e"),
                                   name=NULL, breaks=c("SNP","indel","CDS","UTR"))
  size_scale <- NULL
} else {
  var_layers <- list(
    geom_segment(aes(x=start, xend=xe, y=y, yend=y, colour=vtype), data=vd, linewidth=3.4),
    geom_point(aes(x=xm, y=y, colour=vtype), data=vp, shape=15, size=1.7),
    geom_point(aes(x=start, y=yo, colour=vtype, size=size), data=vi, shape=17))
  size_scale <- scale_size_continuous(range=c(1.8,4), breaks=c(1,5,20), name="insertion (bp)")
  col_scale <- scale_colour_manual(values=c(snp="#2c6fbb", mnp="#7b3294", ins="#1a9850", del="#d6604d",
                                   CDS="#2b2b2b", UTR="#9e9e9e"), name=NULL,
                                   breaks=c("snp","mnp","ins","del","CDS","UTR"),
                                   labels=c("SNP","MNP (2-bp sub)","insertion","deletion","CDS","UTR"))
}

gbot <- p +
  annotate("rect", xmin=HAP_LO, xmax=HAP_HI, ymin=ymin-0.6, ymax=ymax+0.6, fill="#cdd1d5", alpha=0.5) +
  dl() +
  annotate("rect", xmin=xleft, xmax=RHI, ymin=ref_y-0.46, ymax=ref_y+0.46, fill="#cfe8cf", alpha=0.45) +
  annotate("text", x=ref_xc, y=ref_y+0.72, hjust=0.5, size=2.1, fontface="italic", colour="grey15",
           label="reference (TAIR10)") +
  geom_seq(linewidth=0.3, colour="grey82") +
  geom_segment(aes(x=x0, xend=x1, y=gy, yend=gy), data=ext, colour="grey60", linewidth=0.4) +
  geom_segment(aes(x=start, xend=end, y=gy, yend=gy, colour=type), data=filter(mm,type=="UTR"), linewidth=3, lineend="butt") +
  geom_segment(aes(x=start, xend=end, y=gy, yend=gy, colour=type), data=filter(mm,type=="CDS"), linewidth=6, lineend="butt") +
  annotate("text", x=RLO+15, y=gyv+0.6, hjust=0, size=2.6, fontface="italic", colour="grey30", label="AT2G27030.3") +
  annotate("text", x=(11534077+11534176)/2-WIN_LO, y=gyv+0.9, hjust=0.5, size=2.3, fontface="italic", colour="grey30", label="exon 3") +
  var_layers + size_scale + col_scale +
  geom_rect(aes(xmin=xb0, xmax=xb1, ymin=y-0.4, ymax=y+0.4, fill=mean_bio1), data=hap, inherit.aes=FALSE) +
  scale_fill_gradient2(low="#2166ac", mid="#b3b3b3", high="#b2182b", midpoint=B1MEAN, na.value="grey75",
                       name="mean home\nbio1 (C)") +
  geom_text(aes(x=xb0, y=y, label=n_eco), data=hap, hjust=1.2, size=2.5, colour="grey20", inherit.aes=FALSE) +
  geom_text(aes(x=xm, y=glab_top, label=sprintf("%d bp del", size)), data=dlab, size=2.6, colour="#b2182b", fontface="bold") +
  scale_x_continuous(breaks=brk, labels=sprintf("%d", WIN_LO+brk)) +
  coord_cartesian(xlim=c(xleft, RHI), ylim=c(ymin-1, glab_top+1), clip="off") +
  labs(x="Chr2 position (bp)") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

## ---------------- top: LFMM Manhattan, zoomed to the 3' region, frameless ----------------
mv <- av |> mutate(x=pos-WIN_LO, y=-log10(pmax(lfmm_p,1e-300)),
                   shp=ifelse(cls=="snp","SNP","indel")) |> filter(x>=RLO, x<=RHI)
dpmax <- max(abs(mv$dp), na.rm=TRUE); bline <- -log10(bonf)
exon3 <- tibble(s=11534077-WIN_LO, e=11534333-WIN_LO)

gtop <- ggplot(mv, aes(x,y)) +
  annotate("rect", xmin=HAP_LO, xmax=HAP_HI, ymin=-Inf, ymax=Inf, fill="#cdd1d5", alpha=0.5) +
  dl() +
  geom_hline(yintercept=bline, linetype="dashed", colour="grey50") +
  annotate("text", x=RHI, y=bline, label=sprintf("gene-level Bonferroni (0.05/%d)", nrow(av)), hjust=1, vjust=-0.4, size=2.6, colour="grey40") +
  geom_point(aes(colour=dp, shape=shp), size=3.5) +
  scale_colour_gradient2(low="#2166ac", mid="#f7f7f7", high="#b2182b", midpoint=0,
                       limits=c(-dpmax,dpmax), oob=squish, name=expression(Delta*p~"(final - p0)")) +
  scale_shape_manual(values=c(SNP=16, indel=15), name="variant") +
  coord_cartesian(xlim=c(xleft, RHI), clip="off") +
  labs(y=expression(-log[10](p))) +
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.line.x=element_blank())

## ---------------- founder-LD triangle (r^2), hangs below the position axis ----------------
ld <- read_tsv(file.path(HERE,"gg_ld_region.tsv"), show_col_types=FALSE) |> filter(i<=j)
nld <- max(ld$j)+1
gmin <- min(ld$pos_i)-WIN_LO; gmax <- max(ld$pos_j)-WIN_LO; sx <- (gmax-gmin)/(nld-1)
ld <- ld |> mutate(xc=gmin + ((i+j)/2)/(nld-1)*(gmax-gmin), yc=-(j-i)*sx/2)
conn <- data.frame(x=sigx, xend=sigx, y=sx*9, yend=0)   # 3 vertical lines straight down from each dropline
gld <- ggplot(ld, aes(xc, yc, fill=r2)) +
  geom_segment(data=conn, aes(x=x, xend=xend, y=y, yend=yend), linetype="dashed",
               colour="grey45", linewidth=0.5, inherit.aes=FALSE) +
  geom_tile(width=sx*1.04, height=sx*1.04) +
  scale_fill_gradientn(colours=c("#ffffff","#e0e0f0","#b3a2c8","#8073ac","#542788","#2d004b"),
                       limits=c(0,1), name=expression(r^2)) +
  annotate("text", x=barleft, y=0, hjust=0, vjust=1.4, size=2.8, colour="grey35",
           label="founder LD") +
  coord_cartesian(xlim=c(xleft, RHI), ylim=c(min(ld$yc), 0), clip="off") +
  theme_void()

combo <- gtop / gbot / gld + plot_layout(heights=c(1, 1.9, 1.0), guides="collect")
ggsave(file.path(HERE,sprintf("cam5_combined_region%s.png",SUF)), combo, width=12.5, height=12, dpi=170, limitsize=FALSE)
ggsave(file.path(HERE,sprintf("cam5_combined_region%s.pdf",SUF)), combo, width=12.5, height=12, limitsize=FALSE)
cat(sprintf("wrote cam5_combined_region%s.png/.pdf | squares=%s | sig: %s | %d haps\n",
            SUF, SQUARES, paste(sigx, collapse=","), nrow(hap)))
