
tab<- totutable %>%as.data.frame() %>%  rownames_to_column(var="#SampleID") %>% inner_join(metadata)

ro_wet<- tab %>% filter(Type=="Roots", Treatment=="Wet") %>% rename(sampleids="#SampleID") %>% dplyr::select(-BarcodeSequence:-Type)
ro_dry<- tab %>% filter(Type=="Roots", Treatment=="Dry")%>% rename(sampleids="#SampleID")%>% dplyr::select(-BarcodeSequence:-Type)
ro_exdry<- tab %>% filter(Type=="Roots", Treatment=="Extreme_dry")%>% rename(sampleids="#SampleID")%>% dplyr::select(-BarcodeSequence:-Type)

nr_wet<- tab %>% filter(Type=="Bulk soil", Treatment=="Wet")%>% rename(sampleids="#SampleID")%>% dplyr::select(-BarcodeSequence:-Type)
nr_dry<- tab %>% filter(Type=="Bulk soil", Treatment=="Dry")%>% rename(sampleids="#SampleID")%>% dplyr::select(-BarcodeSequence:-Type)
nr_exdry<- tab %>% filter(Type=="Bulk soil", Treatment=="Extreme_dry")%>% rename(sampleids="#SampleID")%>% dplyr::select(-BarcodeSequence:-Type)


ri_wet<- tab %>% filter(Type=="Rhizosphere", Treatment=="Wet")%>% rename(sampleids="#SampleID")%>% dplyr::select(-BarcodeSequence:-Type)
ri_dry<- tab %>% filter(Type=="Rhizosphere", Treatment=="Dry")%>% rename(sampleids="#SampleID")%>% dplyr::select(-BarcodeSequence:-Type)
ri_exdry<- tab %>% filter(Type=="Rhizosphere", Treatment=="Extreme_dry")%>% rename(sampleids="#SampleID")%>% dplyr::select(-BarcodeSequence:-Type)

write_tsv(ro_wet, "table_wet_ro.txt")
write_tsv(ro_dry, "table_dry_ro.txt")
write_tsv(ro_exdry, "table_exdry_ro.txt")

write_tsv(nr_wet, "table_wet_nr.txt")
write_tsv(nr_dry, "table_dry_nr.txt")
write_tsv(nr_exdry, "table_exdry_nr.txt")

write_tsv(ri_wet, "table_wet_ri.txt")
write_tsv(ri_dry, "table_dry_ri.txt")
write_tsv(ri_exdry, "table_exdry_ri.txt")
