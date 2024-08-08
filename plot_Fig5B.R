# Siwei 04 Mar 2024
# Plot Fig.5B

# init ####
{

  library(readxl)
  # library(edgeR)
  library(stringr)
  library(ggplot2)

  # library(data.table)
  library(reshape2)

  library(RColorBrewer)
}

# func ####
# GET EQUATION AND R-SQUARED AS STRING
# SOURCE: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA

lm_eqn <- function(df, x = x, y = y) {
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

# bar_plot #####
df_raw <-
  read_excel(path = "Fig_5B.xlsx")

df_2_plot <-
  df_raw
colnames(df_2_plot) <-
  c("Gene", "qPCR", "RNA-seq")
df_2_plot <-
  melt(df_2_plot,
       value.name = "Gene")
colnames(df_2_plot)[2] <-
  "Measurement"
colnames(df_2_plot)[3] <-
  "value"
df_2_plot$Gene <-
  factor(df_2_plot$Gene,
         levels = unique(df_2_plot$Gene))
## bar plot ####
ggplot(df_2_plot,
       aes(x = Gene,
           y = value,
           group = Measurement,
           fill = Measurement)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = brewer.pal(n = 3,
                                        name = "Set1")[c(2,1)]) +
  ylab("Fold Change") +
  ylim(0, 2.5) +
  # ylim(0, max(df_2_plot$value)) +
  theme_classic() +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"),
        axis.text.x = element_text(angle = 315,
                                   hjust = 0,
                                   vjust = 0))

## correlation ####
df_2_plot <-
  df_raw
colnames(df_2_plot) <-
  c("Gene", "qPCR", "RNA-seq")

ggplot(df_2_plot,
       aes(x = qPCR,
           y = `RNA-seq`)) +
  geom_point() +
  geom_smooth(method = "lm",
              linetype = 2,
              se = F,
              colour = "darkred") +
  scale_x_continuous(limits = c(0, 2.5),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 2.5),
                     expand = c(0, 0)) +
  # xlim(0, 2.5) +
  # ylim(0, 2.5) +
  theme_classic() +
  geom_text(x = 1.7,
            y = 1.8,
            label = lm_eqn(df_2_plot,
                           x = df_2_plot$qPCR,
                           y = df_2_plot$`RNA-seq`),
            parse = T,
            colour = "darkblue") +
  theme(axis.text = element_text(size = 10,
                                 colour = "black"))
