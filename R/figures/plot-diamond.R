# Plot recovery of withheld LC genes using DIAMOnD.
setwd("~/git/spinal-cord-injury-elife-2018")
library(tidyverse)
library(broom)
source("R/theme.R")

# load results
## not commited to git d/t size
results = read.delim("data/PPI/diamond/diamond_final_result.txt.gz")

# summarise bootstraps
sum = results %>%
  group_by(rank, net, type) %>%
  summarise(sampled = mean(sci_sum), random = mean(random_sum), 
            sampled_sd = sd(sci_sum), random_sd = sd(random_sum)) %>%
  gather(recovery, mean, -rank, -net, -type, -sampled_sd, -random_sd) %>%
  filter(type == recovery) %>%
  select(-recovery) %>%
  gather(recovery, sd, -rank, -net, -type, -mean) %>%
  mutate(recovery = gsub("_sd", "", recovery)) %>%
  filter(type == recovery) %>%
  select(-recovery)

# plot Menche interactome for main text
dat = filter(sum, grepl("Menche", net))
p = ggplot(dat, aes(x = rank, y = mean, colour = type)) +
  geom_linerange(color = flat.ui[3], aes(ymin = mean - sd, ymax = mean + sd)) +
  geom_line() +
  labs(y = "Genes recovered", x = "Iterations") +
  scale_color_manual(values = ryb8[c(1, 8)], labels = c("Random", "SCI"), 
                     name = "Seeds") +
  scale_x_continuous(expand = c(0.05, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  sci_theme +
  theme(legend.title = element_text(face = "bold"),
  strip.background = element_rect(fill = 'grey40', color = 'grey40'),
  strip.text = element_text(color = 'white', size = 8),
  axis.line = element_line(color = "grey50"))
p
ggsave("figures/figure-2F.pdf", p, width = 6, height = 4.5, units = 'cm',
       useDingbats = F)

# plot HINT binary interactome
dat = filter(sum, grepl("HINT-b", net))
p = ggplot(dat, aes(x = rank, y = mean, colour = type)) +
  geom_linerange(color = flat.ui[3], aes(ymin = mean - sd, ymax = mean + sd)) +
  geom_line() +
  labs(y = "Genes recovered", x = "Iterations") +
  scale_color_manual(values = ryb8[c(1, 8)], labels = c("Random", "SCI"), 
                     name = "Seeds") +
  scale_x_continuous(expand = c(0.05, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  sci_theme +
  theme(legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = 'grey40', color = 'grey40'),
        strip.text = element_text(color = 'white', size = 8),
        axis.line = element_line(color = "grey50"))
p
ggsave("figures/figure-S1-DIAMOnD-HINT-binary.pdf", p, 
       width = 4.0, height = 4, units = 'cm', useDingbats = F)

# plot HINT cocomplex interactome
dat = filter(sum, grepl("HINT-c", net))
p = ggplot(dat, aes(x = rank, y = mean, colour = type)) +
  geom_linerange(color = flat.ui[3], aes(ymin = mean - sd, ymax = mean + sd)) +
  geom_line() +
  labs(y = "Genes recovered", x = "Iterations") +
  scale_color_manual(values = ryb8[c(1, 8)], labels = c("Random", "SCI"), 
                     name = "Seeds") +
  scale_x_continuous(expand = c(0.05, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  sci_theme +
  theme(legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = 'grey40', color = 'grey40'),
        strip.text = element_text(color = 'white', size = 8),
        axis.line = element_line(color = "grey50"))
p
ggsave("figures/figure-S1-DIAMOnD-HINT-cocomplex.pdf", p, 
       width = 4.0, height = 4, units = 'cm', useDingbats = F)

# plot InBioMap
dat = filter(sum, grepl("InBio", net))
p = ggplot(dat, aes(x = rank, y = mean, colour = type)) +
  geom_linerange(color = flat.ui[3], aes(ymin = mean - sd, ymax = mean + sd)) +
  geom_line() +
  labs(y = "Genes recovered", x = "Iterations") +
  scale_color_manual(values = ryb8[c(1, 8)], labels = c("Random", "SCI"), 
                     name = "Seeds") +
  scale_x_continuous(expand = c(0.05, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  sci_theme +
  theme(legend.title = element_text(face = "bold"),
        strip.background = element_rect(fill = 'grey40', color = 'grey40'),
        strip.text = element_text(color = 'white', size = 8),
        axis.line = element_line(color = "grey50"))
p
ggsave("figures/figure-S1-DIAMOnD-InBioMap.pdf", p, 
       width = 4.0, height = 4, units = 'cm', useDingbats = F)

# do K-S tests
sum %>%
  select(-sd) %>%
  group_by(rank, net) %>%
  spread(type, mean) %>%
  group_by(net) %>%
  do(ks = ks.test(.$sampled, .$random)) %>%
  tidy(ks)
