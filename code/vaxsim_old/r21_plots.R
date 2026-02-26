fig_stochastic_analytic_m72IIb_r21 <- plot_stochastic_analytic(stochastic_df_m72IIb, analytical_df_m72IIb, cols=c("n_tested","n_recruited")) + 
  labs(title="M72/AS01E trial emulation", subtitle=paste0("age_structure = ", age_structure))
fig_stochastic_analytic_m72IIb_r21
ggsave(fig_stochastic_analytic_m72IIb_r21, file=paste0("figures/stochastic_analytic_m72IIb_r21_",age_structure,".pdf"), width=5, height=5/1.6)

# fig_stochastic_analytic_m72IIb_r21_log <- plot_stochastic_analytic(stochastic_df_m72IIb, analytical_df_m72IIb, 
#                                                    cols=c("n_tested","n_recruited")) + 
# 	labs(title="M72/AS01E trial emulation") + 
# 	scale_y_continuous(trans="log10")
# fig_stochastic_analytic_m72IIb_r21_log
# ggsave(fig_stochastic_analytic_m72IIb_r21_log, file="figures/stochastic_analytic_m72IIb_r21_log.pdf", width=5, height=5/1.6)