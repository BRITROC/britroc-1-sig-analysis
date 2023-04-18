patient.meta[patient.meta$PATIENT_ID == 'BRITROC-67',]

samples_interest <- patient.meta$SAMPLE_ID[(patient.meta$PATIENT_ID %in% c('BRITROC-23', 'BRITROC-74', 'BRITROC-209', 'BRITROC-216', 
'BRITROC-241', 'BRITROC-267', 'BRITROC-274')) & (patient.meta$group == 'arx') ]

exp_arx_undergoingWGD <- exposures[as.character(samples_interest),]
exp_arx_undergoingWGD
createBarplot(exp_arx_undergoingWGD)

patient.meta[patient.meta$PATIENT_ID == 'BRITROC-274',]
patient.meta[patient.meta$PATIENT_ID == 'BRITROC-232',]

patient.meta$PATIENT_ID <- as.character(patient.meta$PATIENT_ID)
patient.meta$SAMPLE_ID <- as.character(patient.meta$SAMPLE_ID)
patient.meta$SAMPLE_ID_long <- paste0(patient.meta$SAMPLE_ID, '_', patient.meta$group)

multiple_samples_at_least_one <- patient.meta %>% group_by(PATIENT_ID) %>% summarise(tab=length(unique(group))) %>% filter(tab == 2)
multiple_samples_at_least_one <- patient.meta[match(unlist(multiple_samples_at_least_one[,1]), patient.meta$PATIENT_ID),]
length(unique(multiple_samples_at_least_one$PATIENT_ID))

multiple_samples_at_least_one[multiple_samples_at_least_one$group == "arx",]

multiple_samples <- patient.meta %>% group_by(PATIENT_ID) %>% summarise(tab=length(unique(group)), all=length(group)) %>% 
  filter(tab == 2, all>2)
multiple_samples

rownames(exposures) <- patient.meta$SAMPLE_ID_long[match(rownames(exposures), patient.meta$SAMPLE_ID)]

remove_sig <- function(exp, sig_to_remove){
  normalise_rw(exp[,! (colnames(exp) %in% sig_to_remove)])
}

patient.meta[patient.meta$PATIENT_ID == as.character(multiple_samples[1,1]),]


pdf("../../../results/exploratory/multiple_samples_exposures.pdf")
for( mult_patient in 1:nrow(multiple_samples)){
  patient_name <- as.character(multiple_samples[mult_patient,1])
  print(createBarplot(exposures[patient.meta$SAMPLE_ID_long[which(patient.meta$PATIENT_ID == patient_name)],],
                angle_rotation_axis = 45)+ggtitle(patient_name))
}
dev.off()


# createBarplot(remove_sig(exposures[patient.meta$SAMPLE_ID_long[which(patient.meta$PATIENT_ID == as.character(multiple_samples[1,1]))],],
#               's5'))
# createBarplot(remove_sig(exposures[patient.meta$SAMPLE_ID_long[which(patient.meta$PATIENT_ID == as.character(multiple_samples[1,1]))],],
#                          's1'))
