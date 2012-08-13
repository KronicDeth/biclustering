load PreGeneGroups;
load PretrimList14;

Biogen=sortrows(([Biogen;Unclass;321]));
BiogenSigUn=significant(Biogen,Genes_in_Patterng);
save PreBioList14TX BiogenSigUn 

CellCom=sortrows(([CellCom;Unclass;321]));
CellComSigUn=significant(CellCom,Genes_in_Patterng);
save PreBioList14TX CellComSigUn -append

CellCycle=sortrows(([CellCycle;Unclass;321]));
CellCycleSigUn=significant(CellCycle,Genes_in_Patterng);
save PreBioList14TX CellCycleSigUn -append

CellEnv=sortrows(([CellEnv;Unclass;321]));
CellEnvSigUn=significant(CellEnv,Genes_in_Patterng);
save PreBioList14TX CellEnvSigUn -append

CellRescue=sortrows(([CellRescue;Unclass;321]));
CellRescueSigUn=significant(CellRescue,Genes_in_Patterng);
save PreBioList14TX CellRescueSigUn -append

CellTrans=sortrows(([CellTrans;Unclass;321]));
CellTransSigUn=significant(CellTrans,Genes_in_Patterng);
save PreBioList14TX CellTransSigUn -append

CellType=sortrows(([CellType;Unclass;321]));
CellTypeSigUn=significant(CellType,Genes_in_Patterng);
save PreBioList14TX CellTypeSigUn -append

Energy=sortrows(([Energy;Unclass;321]));
EnergySigUn=significant(Energy, Genes_in_Patterng);
save PreBioList14TX EnergySigUn -append

Metabolism=sortrows(([Metabolism;Unclass;321]));
MetabolismSigUn=significant(Metabolism,Genes_in_Patterng);
save PreBioList14TX MetabolismSigUn -append

ProteinAct=sortrows(([ProteinAct;Unclass;321]));
ProteinActSigUn=significant(ProteinAct,Genes_in_Patterng);
save PreBioList14TX ProteinActSigUn -append

ProteinBind=sortrows(([ProteinBind;Unclass;321]));
ProteinBindSigUn=significant(ProteinBind, Genes_in_Patterng);
save PreBioList14TX ProteinBindSigUn -append

ProteinFate=sortrows(([ProteinFate;Unclass;321]));
ProteinFateSigUn=significant(ProteinFate,Genes_in_Patterng);
save PreBioList14TX ProteinFateSigUn -append

ProteinSynth=sortrows(([ProteinSynth;Unclass;321]));
ProteinSynthSigUn=significant(ProteinSynth,Genes_in_Patterng);
save PreBioList14TX ProteinSynthSigUn -append

SysDev=sortrows(([SysDev;Unclass;321]));
SysDevSigUn=significant(SysDev,Genes_in_Patterng);
save PreBioList14TX SysDevSigUn -append

SysEnv=sortrows(([SysDev;Unclass;321]));
SysEnvSigUn=significant(SysEnv,Genes_in_Patterng);
save PreBioList14TX SysEnvSigUn -append

Transcription=sortrows(([Transcription;Unclass;321]));
TranscriptionSigUn=significant(Transcription,Genes_in_Patterng);
save PreBioList14TX TranscriptionSigUn -append

Transposable=sortrows(([Transposable;Unclass;321]));
TransposableSigUn=significant(Transposable,Genes_in_Patterng);
save PreBioList14TX TransposableSigUn -append