ggplot(data = cdp, aes(x = CONDITION, y = WTKg)) + 
  geom_boxplot(lwd = 0.5, na.rm = TRUE) + 
  geom_point(size = 2, position = position_dodge2(0.1), na.rm = TRUE) + 
  stat_compare_means(comparisons = comparison1, size = 3, na.rm = FALSE,
                     method = "t.test", paired = TRUE, label = "p.format") + 
  stat_compare_means(comparisons = comparison2, size = 3, na.rm = FALSE,
                     method = "t.test", paired = FALSE, label = "p.format") +
  ylab("Weight(Kg)") +
  xlab("") + 
  theme_bw() + 
  theme(text = element_text(size = 15)) + 
  ggsave(filename = "fig1.png", device = "png", dpi = 1200, 
         width = 6, height = 5)

ggplot(data = cdp, aes(x = CONDITION, y = Hb)) + 
  geom_boxplot(lwd = 0.5, na.rm = TRUE) + 
  geom_point(size = 2, position = position_dodge2(0.1), na.rm = TRUE) + 
  stat_compare_means(comparisons = comparison1, size = 3, na.rm = FALSE,
                     method = "t.test", paired = TRUE, label = "p.format") + 
  stat_compare_means(comparisons = comparison2, size = 3, na.rm = FALSE,
                     method = "t.test", paired = FALSE, label = "p.format") +
  ylab("[Hb] (g/dL)") +
  xlab("") + 
  theme_bw() + 
  theme(text = element_text(size = 15)) + 
  ggsave(filename = "fig2.png", device = "png", dpi = 1200, 
         width = 6, height = 5)

ggplot(data = cdp, aes(x = CONDITION, y = P50)) + 
  geom_boxplot(lwd = 0.5, na.rm = TRUE) + 
  geom_point(size = 2, position = position_dodge2(0.1), na.rm = TRUE) + 
  stat_compare_means(comparisons = comparison1, size = 3, na.rm = FALSE,
                     method = "t.test", paired = TRUE, label = "p.format") + 
  stat_compare_means(comparisons = comparison2, size = 3, na.rm = FALSE,
                     method = "t.test", paired = FALSE, label = "p.format") +
  ylab(TeX("$P_{50}$")) +
  xlab("") + 
  theme_bw() + 
  theme(text = element_text(size = 15)) + 
  ggsave(filename = "fig3.png", device = "png", dpi = 1200, 
         width = 6, height = 5)

ggplot(data = cdp, aes(x = CONDITION, y = CaO2)) + 
  geom_boxplot(lwd = 0.5, na.rm = TRUE) + 
  geom_point(size = 2, position = position_dodge2(0.1), na.rm = TRUE) + 
  stat_compare_means(comparisons = comparison1, size = 3, na.rm = FALSE,
                     method = "t.test", paired = TRUE, label = "p.format") + 
  stat_compare_means(comparisons = comparison2, size = 3, na.rm = FALSE,
                     method = "t.test", paired = FALSE, label = "p.format") +
  ylab(TeX("$C_{a}O_{2}$")) +
  xlab("") + 
  theme_bw() + 
  theme(text = element_text(size = 15)) + 
  ggsave(filename = "fig4.png", device = "png", dpi = 1200, 
         width = 6, height = 5)

ggplot(data = cdp, aes(x = CONDITION, y = CVO2)) + 
  geom_boxplot(lwd = 0.5, na.rm = TRUE) + 
  geom_point(size = 2, position = position_dodge2(0.1), na.rm = TRUE) + 
  stat_compare_means(comparisons = comparison1, size = 3, na.rm = FALSE,
                     method = "t.test", paired = TRUE, label = "p.format") + 
  stat_compare_means(comparisons = comparison2, size = 3, na.rm = FALSE,
                     method = "t.test", paired = FALSE, label = "p.format") +
  ylab(TeX("$C_{V}O_{2}$")) +
  xlab("") + 
  theme_bw() + 
  theme(text = element_text(size = 15)) + 
  ggsave(filename = "fig5.png", device = "png", dpi = 1200, 
         width = 6, height = 5)

ggplot(data = cdp, aes(x = CONDITION, y = VO2maxPerKg)) + 
  geom_boxplot(lwd = 0.5, na.rm = TRUE) + 
  geom_point(size = 2, position = position_dodge2(0.1), na.rm = TRUE) + 
  stat_compare_means(comparisons = comparison1, size = 3, na.rm = FALSE,
                     method = "t.test", paired = TRUE, label = "p.format") + 
  stat_compare_means(comparisons = comparison2, size = 3, na.rm = FALSE,
                     method = "t.test", paired = FALSE, label = "p.format") +
  ylab(TeX("$VO_{2max} Kg^{-1}$")) +
  xlab("") + 
  theme_bw() + 
  theme(text = element_text(size = 15)) + 
  ggsave(filename = "fig6.png", device = "png", dpi = 1200, 
         width = 6, height = 5)

ggplot(data = cdp, aes(x = CONDITION, y = QTmaxPerKg)) + 
  geom_boxplot(lwd = 0.5, na.rm = TRUE) + 
  geom_point(size = 2, position = position_dodge2(0.1), na.rm = TRUE) + 
  stat_compare_means(comparisons = comparison1, size = 3, na.rm = FALSE,
                     method = "t.test", paired = TRUE, label = "p.format") + 
  stat_compare_means(comparisons = comparison2, size = 3, na.rm = FALSE,
                     method = "t.test", paired = FALSE, label = "p.format") +
  ylab(TeX("$QT_{max} Kg^{-1}$")) +
  xlab("") + 
  theme_bw() + 
  theme(text = element_text(size = 15)) + 
  ggsave(filename = "fig7.png", device = "png", dpi = 1200, 
         width = 6, height = 5)
