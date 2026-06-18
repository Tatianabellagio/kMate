#!/usr/bin/env Rscript
# Combined CAM5 figure: LFMM GEA Manhattan (top) + unique-haplotype pangenome
# panel (bottom), sharing the genomic x-axis. Dashed droplines connect the
# Bonferroni-significant LFMM variants down into the haplotype map.
# Arg "squares": draw variants as SNP=circle / indel=square (matching the
# Manhattan) instead of type/size glyphs.
suppressMessages({library(gggenomes); library(dplyr); library(readr)
                  library(patchwork); library(scales)})
HERE <- "/global/scratch/users/tbellg/kmate/analysis/grenenet_gea/cam5_replication/pang_cam5"
infct <- function(x) factor(x, levels = unique(x))
LEN <- 2533; WIN_LO <- 11531800
HAP_LO <- 11533904-WIN_LO; HAP_HI <- 11534143-WIN_LO
SQUARES <- length(commandArgs(TRUE))>0 && commandArgs(TRUE)[1]=="squares"
SUF <- if (SQUARES) "_squares" else ""

# significant LFMM variants + dropline helper
av <- read_tsv(file.path(HERE,"cam5_variants_gen9.tsv"), show_col_types=FALSE) |> filter(!is.na(lfmm_p))
bonf <- 0.05/nrow(av)
sigx <- (av |> filter(lfmm_p<bonf) |> mutate(x=pos-WIN_LO))$x
dl <- function() geom_vline(xintercept=sigx, linetype="dashed", colour="grey45", linewidth=0.45, alpha=0.55)

## ---------------- bottom panel: haplotype pangenome ----------------
seqs <- read_tsv(file.path(HERE,"gg_hapseqs.tsv"),  show_col_types=FALSE) |>
  mutate(seq_id=as.character(seq_id)) |> filter(kind=="hap") |> mutate(bin_id=infct(seq_id))
model <- read_tsv(file.path(HERE,"gg_model.tsv"),  show_col_types=FALSE) |> mutate(seq_id=as.character(seq_id))
vars  <- read_tsv(file.path(HERE,"gg_hapfeats.tsv"), show_col_types=FALSE) |> mutate(seq_id=as.character(seq_id))
B1MEAN <- read_tsv(file.path(HERE,"gg_hapmeta.tsv"), show_col_types=FALSE)$x[1]

p  <- gggenomes(seqs=seqs); ly <- pull_seqs(p); hap <- ly
ymin <- min(ly$y); ymax <- max(ly$y)
v  <- vars |> left_join(select(ly, seq_id, y), by="seq_id") |> filter(!is.na(y))
ilab <- v |> filter(size>=20) |> distinct(start, end, vtype, size) |>
  mutate(xm=(start+end)/2, lab=sprintf("%d bp %s", size, ifelse(vtype=="del","del","ins")))
ISO <- "AT2G27030.3"; gyv <- ymax + 2.6
ext <- model |> filter(seq_id==ISO) |> summarise(x0=min(start), x1=max(end)) |> mutate(gy=gyv)
mm  <- model |> filter(seq_id==ISO) |> mutate(gy=gyv)
glab_top <- ymax + 4.2
BARGAP <- 10; BARMAX <- 95; K <- BARMAX/sqrt(max(hap$n_eco))
hap <- hap |> mutate(xb1=-BARGAP, xb0=-BARGAP - sqrt(n_eco)*K)
barleft <- -BARGAP - BARMAX - 20; xleft <- barleft
ref_id <- setdiff(seqs$seq_id, unique(vars$seq_id))[1]
ref_y  <- ly$y[ly$seq_id==ref_id]; ref_n <- hap$n_eco[hap$seq_id==ref_id]
hap <- hap |> mutate(lab = ifelse(seq_id==ref_id, sprintf("reference  %d", n_eco), as.character(n_eco)))
brk <- c(0,500,1000,1500,2000,2500)

# variant layers + colour scale depend on mode
if (SQUARES) {
  vv <- v |> mutate(grp=ifelse(vtype=="snp","SNP","indel"), xm=(start+end)/2)
  var_layers <- list(geom_point(aes(x=xm, y=y, colour=grp, shape=grp), data=vv, size=1.6),
                     scale_shape_manual(values=c(SNP=16, indel=15), name="variant"))
  col_scale <- scale_colour_manual(values=c(SNP="grey30", indel="grey55", CDS="#2b2b2b", UTR="#9e9e9e"),
                                   name=NULL, breaks=c("SNP","indel","CDS","UTR"))
  size_scale <- NULL
} else {
  vp <- v |> filter(vtype %in% c("snp","mnp")) |> mutate(xm=(start+end)/2)
  vi <- v |> filter(vtype=="ins") |> mutate(yo=y+0.30)
  vd <- v |> filter(vtype=="del") |> mutate(xe=pmax(end, start+2))
  var_layers <- list(
    geom_segment(aes(x=start, xend=xe, y=y, yend=y, colour=vtype), data=vd, linewidth=2.7),
    geom_point(aes(x=xm, y=y, colour=vtype), data=vp, shape=15, size=1.25),
    geom_point(aes(x=start, y=yo, colour=vtype, size=size), data=vi, shape=17))
  size_scale <- scale_size_continuous(range=c(1.5,4), breaks=c(1,50,500,2668), name="insertion (bp)")
  col_scale <- scale_colour_manual(values=c(snp="#2c6fbb", mnp="#7b3294", ins="#1a9850", del="#d6604d",
                                   CDS="#2b2b2b", UTR="#9e9e9e"), name=NULL,
                                   breaks=c("snp","mnp","ins","del","CDS","UTR"),
                                   labels=c("SNP","MNP (2-bp sub)","insertion","deletion","CDS","UTR"))
}

gbot <- p +
  annotate("rect", xmin=HAP_LO, xmax=HAP_HI, ymin=ymin-0.6, ymax=ymax+0.6, fill="#cdd1d5", alpha=0.5) +
  dl() +
  annotate("rect", xmin=xleft, xmax=LEN, ymin=ref_y-0.46, ymax=ref_y+0.46, fill="#cfe8cf", alpha=0.5) +
  geom_seq(linewidth=0.18, colour="grey85") +
  geom_segment(aes(x=x0, xend=x1, y=gy, yend=gy), data=ext, colour="grey60", linewidth=0.4) +
  geom_segment(aes(x=start, xend=end, y=gy, yend=gy, colour=type), data=filter(mm,type=="UTR"), linewidth=2.6, lineend="butt") +
  geom_segment(aes(x=start, xend=end, y=gy, yend=gy, colour=type), data=filter(mm,type=="CDS"), linewidth=5, lineend="butt") +
  geom_text(aes(x=x0, y=gy+1.0, label="AT2G27030.3"), data=ext, hjust=0, size=2.6, colour="grey30", fontface="italic") +
  var_layers + size_scale + col_scale +
  geom_rect(aes(xmin=xb0, xmax=xb1, ymin=y-0.42, ymax=y+0.42, fill=mean_bio1), data=hap, inherit.aes=FALSE) +
  scale_fill_gradient2(low="#2166ac", mid="#b3b3b3", high="#b2182b", midpoint=B1MEAN, na.value="grey75",
                       name="mean home\nbio1 (C)") +
  geom_text(aes(x=xb0, y=y, label=lab), data=filter(hap, n_eco>=4 | seq_id==ref_id), hjust=1.15, size=1.9, colour="grey25", inherit.aes=FALSE) +
  geom_text(aes(x=xm, y=glab_top, label=lab), data=ilab, size=2.3, colour="#8c2d04", fontface="bold") +
  scale_x_continuous(breaks=brk, labels=sprintf("%d", WIN_LO+brk)) +
  coord_cartesian(xlim=c(xleft, LEN), ylim=c(ymin-1, glab_top+1), clip="off") +
  labs(x="Chr2 position (bp)") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

## ---------------- top panel: LFMM Manhattan ----------------
mv <- av |> mutate(x=pos-WIN_LO, y=-log10(pmax(lfmm_p,1e-300)), shp=ifelse(cls=="snp","SNP","indel"))
dpmax <- max(abs(mv$dp), na.rm=TRUE); bline <- -log10(bonf)
exons <- tibble(s=c(11532004,11532687,11534077)-WIN_LO, e=c(11532144,11533056,11534333)-WIN_LO)

gtop <- ggplot(mv, aes(x,y)) +
  annotate("rect", xmin=HAP_LO, xmax=HAP_HI, ymin=-Inf, ymax=Inf, fill="#cdd1d5", alpha=0.5) +
  dl() +
  geom_hline(yintercept=bline, linetype="dashed", colour="grey50") +
  annotate("text", x=LEN, y=bline, label=sprintf("gene-level Bonferroni (0.05/%d)", nrow(av)), hjust=1, vjust=-0.4, size=2.6, colour="grey40") +
  geom_point(aes(colour=dp, shape=shp), size=2.8) +
  scale_colour_gradient2(low="#2166ac", mid="#f7f7f7", high="#b2182b", midpoint=0,
                       limits=c(-dpmax,dpmax), oob=squish, name=expression(Delta*p~"(final - p0)")) +
  scale_shape_manual(values=c(SNP=16, indel=15), name="variant") +
  coord_cartesian(xlim=c(xleft, LEN), clip="off") +
  labs(y=expression(-log[10](p))) +
  theme_classic() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.line.x=element_blank())

## ---------------- founder-LD triangle (r^2), hangs below the position axis ----------------
ld <- read_tsv(file.path(HERE,"gg_ld_gene.tsv"), show_col_types=FALSE) |> filter(i<=j)
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
  coord_cartesian(xlim=c(xleft, LEN), ylim=c(min(ld$yc), 0), clip="off") +
  theme_void()

combo <- gtop / gbot / gld + plot_layout(heights=c(1, 3.3, 1.2), guides="collect")
ggsave(file.path(HERE,sprintf("cam5_combined_manhattan_hap%s.png",SUF)), combo, width=13, height=19, dpi=160, limitsize=FALSE)
ggsave(file.path(HERE,sprintf("cam5_combined_manhattan_hap%s.pdf",SUF)), combo, width=13, height=19, limitsize=FALSE)
cat(sprintf("wrote cam5_combined_manhattan_hap%s.png/.pdf | squares=%s | sig=%s\n",
            SUF, SQUARES, paste(sigx, collapse=",")))
