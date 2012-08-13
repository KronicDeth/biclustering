load List16;
load RegBioListT16;

Bio=BiogenSig.*TGenes_in_Patterng;
CCom=CellComSig.*TGenes_in_Patterng;
CCycle=CellCycleSig.*TGenes_in_Patterng;
CEnv=CellEnvSig.*TGenes_in_Patterng;
CRescue=CellRescueSig.*TGenes_in_Patterng;
CTrans=CellTransSig.*TGenes_in_Patterng;
CType=CellTypeSig.*TGenes_in_Patterng;
Energy=EnergySig.*TGenes_in_Patterng;
Meta=MetabolismSig.*TGenes_in_Patterng;
PAct=ProteinActSig.*TGenes_in_Patterng;
PBind=ProteinBindSig.*TGenes_in_Patterng;
PFate=ProteinFateSig.*TGenes_in_Patterng;
PSynth=ProteinSynthSig.*TGenes_in_Patterng;
SDev=SysDevSig.*TGenes_in_Patterng;
SEnv=SysEnvSig.*TGenes_in_Patterng;
Transcript=TranscriptionSig.*TGenes_in_Patterng;
Transpos=TransposableSig.*TGenes_in_Patterng;


save BreakdownReg16 Bio CCom CCycle CEnv CRescue CTrans CType Energy Meta PAct PBind PFate PSynth SDev SEnv Transcript Transpos