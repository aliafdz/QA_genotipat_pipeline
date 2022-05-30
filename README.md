# QA and genotyping pipeline

## Descripció
En aquest repositori es troba el *pipeline* original emprat per a la realització del Treball de Fi de Màster de l'alumna Alicia Aranda Fernandez, estudiant del Màster en Bioinformàtica i Bioestadística (UOC-UB). 

## Estructura del directori
Aquest directori només inclou la part corresponent als fitxers de codi per tal de realitzar un seguiment sobre la seva revisió, que ha consisit en documentar tots els passos realitzats al llarg de cada *script*. 

La carpeta [R](./R) inclou tots els *scripts* que es van cridant al llarg dels arxius globals [MiSeq_QA_Pipeline-v2.5.R](./MiSeq_QA_Pipeline-v2.5.R) i [HBV_Genotype_MySeq_Pipeline-v1.3.R](./HBV_Genotype_MySeq_Pipeline-v1.3.R). Cal tenir en compte, però, que només una secció del *pipeline* proporcionat ha estat revisada, per la qual cosa els arxius que no han estat emprats en la consecució del projecte es troben en la subcarpeta [ignore](./R/ignore). 

La carpeta [VHBass3_pipeline](./VHBass3_pipeline) està constituida per les subcarpetes [reports](./VHBass3_pipeline/reports) i [results](./VHBass3_pipeline/results), i emmagatzema tots els arxius resultants de l'aplicació del *pipeline* original al complet sobre el conjunt de dades identificades amb el nom de projecte `VHBass3` (dades assistencials). Alguns dels gràfics presentats també es descriuen a la memòria del TFM corresponent. 

## Autor/a
El codi original va ser dissenyat pel Dr. Josep Gregori i Font del [grup en Malalties Hepàtiques del VHIR](<http://www.vhir.org/web_vhir/portal1/grup-presentacio.asp?s=recerca&contentid=187009&idrefer=187010>).

Les revisions i modificacions pertinents del codi han estat realitzades per Alicia Aranda Fernandez.