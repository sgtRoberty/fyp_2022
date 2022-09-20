load("../results/mccr_test/malaysia_ltt_mccr_output.RData")
load("../results/mccr_test/panama_ltt_mccr_output.RData")

# Malaysia
malaysia_ltt
malaysia_mccr_1000sim
malaysia_mccr_5000sim
plot(malaysia_mccr_1000sim)
plot(malaysia_mccr_5000sim)

malaysia_mccrTest_1000sim
malaysia_mccrTest_5000sim

mean(malaysia_mccr_1000sim$null.gamma)
mean(malaysia_mccrTest_1000sim$null.gamma)

# Panama
panama_ltt
panama_mccr_1000sim
panama_mccr_5000sim
plot(panama_mccr_1000sim)
plot(panama_mccr_5000sim)

panama_mccrTest_1000sim
panama_mccrTest_5000sim

mean(panama_mccr_1000sim$null.gamma)
mean(panama_mccrTest_1000sim$null.gamma)
