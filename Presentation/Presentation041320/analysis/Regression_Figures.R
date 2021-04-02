ggplot(data = subset(cdp, HMD == 0), aes(x = Hb, y = P50, group = CMS)) + 
  geom_point(aes(color = CONDITION), size = 2) + 
  geom_smooth(method = "lm",formula = "y~x", color = "black", size = 1.2) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = "y~x"),
                  geom = "text",
                  na.rm = TRUE,
                  aes(label = paste("P = ", signif(..p.value.., digits = 3), sep = "")),
                  size = 3) +
  xlab("[Hb] (g/dL)") + 
  ylab(TeX("$P_{50}$")) + 
  theme_bw() +
  theme(text = element_text(size = 15)) +
  ggsave(filename = "fig8.png", device = "png", dpi = 1200, 
         width = 7, height = 5)

ggplot(data = subset(cdp, HMD == 0), aes(x = Hb, y = CaO2, group = CMS)) + 
  geom_point(aes(color = CONDITION), size = 2) + 
  geom_smooth(method = "lm",formula = "y~x", color = "black", size = 1.2) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = "y~x"),
                  geom = "text",
                  na.rm = TRUE,
                  aes(label = paste("P = ", signif(..p.value.., digits = 3), sep = "")),
                  size = 3) +
  xlab("[Hb] (g/dL)") + 
  ylab(TeX("$C_{a}O_{2}$")) + 
  theme_bw() +
  theme(text = element_text(size = 15)) +
  ggsave(filename = "fig9.png", device = "png", dpi = 1200, 
         width = 7, height = 5)

ggplot(data = subset(cdp, HMD == 0), aes(x = QTmaxPerKg, y = VO2maxPerKg, group = CMS)) + 
  geom_point(aes(color = CONDITION), size = 2) + 
  geom_smooth(method = "lm",formula = "y~x", color = "black", size = 1.2) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = "y~x"),
                  geom = "text",
                  na.rm = TRUE,
                  aes(label = paste("P = ", signif(..p.value.., digits = 3), sep = "")),
                  size = 3) +
  xlab(TeX("$QT_{max} Kg^{-1}$")) +
  ylab(TeX("$VO_{2max} Kg^{-1}$")) +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  ggsave(filename = "fig10.png", device = "png", dpi = 1200, 
         width = 7, height = 5)

ggplot(data = subset(cdp, HMD == 0), aes(x = WTKg, y = VO2maxPerKg, group = CMS)) + 
  geom_point(aes(color = CONDITION), size = 2) + 
  geom_smooth(method = "lm",formula = "y~x", color = "black", size = 1.2) +
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = "y~x"),
                  geom = "text",
                  na.rm = TRUE,
                  aes(label = paste("P = ", signif(..p.value.., digits = 3), sep = "")),
                  size = 3) +
  xlab(TeX("$Weight(Kg)$")) +
  ylab(TeX("$VO_{2max} Kg^{-1}$")) +
  theme_bw() +
  theme(text = element_text(size = 15)) +
  ggsave(filename = "fig11.png", device = "png", dpi = 1200, 
         width = 7, height = 5)