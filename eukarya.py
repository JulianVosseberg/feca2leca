#!/usr/bin/python3

# Supergroup variables
supergroups5 = {'ACAS' : 'Amoebozoa',
                'EHIS' : 'Amoebozoa',
                'EINV' : 'Amoebozoa',
                'DDIS' : 'Amoebozoa',
                'PPAL' : 'Amoebozoa',
                'ASUB' : 'Amoebozoa',
                'TTRA' : 'Obazoa',
                'SARC' : 'Obazoa',
                'CFRA' : 'Obazoa',
                'COWC' : 'Obazoa',
                'MBRE' : 'Obazoa',
                'SROS' : 'Obazoa',
                'TADH' : 'Obazoa',
                'AQUE' : 'Obazoa',
                'ADIG' : 'Obazoa',
                'NVEC' : 'Obazoa',
                'HVUL' : 'Obazoa',
                'TKIT' : 'Obazoa',
                'MLEI' : 'Obazoa',
                'AVAG' : 'Obazoa',
                'EMUL' : 'Obazoa',
                'TSOL' : 'Obazoa',
                'TASI' : 'Obazoa',
                'GSAL' : 'Obazoa',
                'SMED' : 'Obazoa',
                'SMAN' : 'Obazoa',
                'PCAU' : 'Obazoa',
                'SRAT' : 'Obazoa',
                'SBAT' : 'Obazoa',
                'OVOL' : 'Obazoa',
                'PPAC' : 'Obazoa',
                'BMAL' : 'Obazoa',
                'CELE' : 'Obazoa',
                'TSPI' : 'Obazoa',
                'BANT' : 'Obazoa',
                'DMEL' : 'Obazoa',
                'AGAM' : 'Obazoa',
                'GMOR' : 'Obazoa',
                'BMOR' : 'Obazoa',
                'APIS' : 'Obazoa',
                'TCAS' : 'Obazoa',
                'RPRO' : 'Obazoa',
                'ZNEV' : 'Obazoa',
                'PHUM' : 'Obazoa',
                'ISCA' : 'Obazoa',
                'SMIM' : 'Obazoa',
                'SMAR' : 'Obazoa',
                'LSAL' : 'Obazoa',
                'DMAG' : 'Obazoa',
                'DPUL' : 'Obazoa',
                'LANA' : 'Obazoa',
                'HAZT' : 'Obazoa',
                'LPOL' : 'Obazoa',
                'RVAR' : 'Obazoa',
                'CGIG' : 'Obazoa',
                'LGIG' : 'Obazoa',
                'BGLA' : 'Obazoa',
                'OBIM' : 'Obazoa',
                'HROB' : 'Obazoa',
                'CTEL' : 'Obazoa',
                'ILIN' : 'Obazoa',
                'PFLA' : 'Obazoa',
                'SKOW' : 'Obazoa',
                'PMIN' : 'Obazoa',
                'SPUR' : 'Obazoa',
                'BBEL' : 'Obazoa',
                'BFLO' : 'Obazoa',
                'ODIO' : 'Obazoa',
                'CINT' : 'Obazoa',
                'BSCH' : 'Obazoa',
                'PETM' : 'Obazoa',
                'CMIL' : 'Obazoa',
                'DRER' : 'Obazoa',
                'TRUB' : 'Obazoa',
                'XTRO' : 'Obazoa',
                'APLA' : 'Obazoa',
                'PSIN' : 'Obazoa',
                'OANA' : 'Obazoa',
                'MMUS' : 'Obazoa',
                'HSAP' : 'Obazoa',
                'NUSP' : 'Obazoa',
                'FALB' : 'Obazoa',
                'RALL' : 'Obazoa',
                'MDAP' : 'Obazoa',
                'VCUL' : 'Obazoa',
                'EINT' : 'Obazoa',
                'PSPF' : 'Obazoa',
                'PSPE' : 'Obazoa',
                'NSPE' : 'Obazoa',
                'ASPE' : 'Obazoa',
                'OSPE' : 'Obazoa',
                'BDEN' : 'Obazoa',
                'SPUN' : 'Obazoa',
                'GPRO' : 'Obazoa',
                'AMAC' : 'Obazoa',
                'CANG' : 'Obazoa',
                'BBRI' : 'Obazoa',
                'BMER' : 'Obazoa',
                'CCOR' : 'Obazoa',
                'CREV' : 'Obazoa',
                'SMUC' : 'Obazoa',
                'LPEN' : 'Obazoa',
                'MPTE' : 'Obazoa',
                'RBRE' : 'Obazoa',
                'SPLU' : 'Obazoa',
                'MCIR' : 'Obazoa',
                'PBLA' : 'Obazoa',
                'LCOR' : 'Obazoa',
                'LTRA' : 'Obazoa',
                'MVER' : 'Obazoa',
                'MELO' : 'Obazoa',
                'RIRR' : 'Obazoa',
                'RSPE' : 'Obazoa',
                'UMAY' : 'Obazoa',
                'AING' : 'Obazoa',
                'MLAR' : 'Obazoa',
                'MOSM' : 'Obazoa',
                'ABIS' : 'Obazoa',
                'RSOL' : 'Obazoa',
                'CNEO' : 'Obazoa',
                'SCOM' : 'Obazoa',
                'SPOM' : 'Obazoa',
                'PMUR' : 'Obazoa',
                'SCER' : 'Obazoa',
                'KLAC' : 'Obazoa',
                'YLIP' : 'Obazoa',
                'NCRA' : 'Obazoa',
                'FOXY' : 'Obazoa',
                'TMEL' : 'Obazoa',
                'BHOM' : 'RASH',
                'AKER' : 'RASH',
                'SAGG' : 'RASH',
                'ALIM' : 'RASH',
                'PINF' : 'RASH',
                'HPAR' : 'RASH',
                'PHAL' : 'RASH',
                'PULT' : 'RASH',
                'ALAI' : 'RASH',
                'AAST' : 'RASH',
                'SPAR' : 'RASH',
                'NGAD' : 'RASH',
                'AANO' : 'RASH',
                'ESIL' : 'RASH',
                'COKA' : 'RASH',
                'PTRI' : 'RASH',
                'TPSE' : 'RASH',
                'FCYL' : 'RASH',
                'PMAR' : 'RASH',
                'SMIC' : 'RASH',
                'SMIN' : 'RASH',
                'PFAL' : 'RASH',
                'BBIG' : 'RASH',
                'TANN' : 'RASH',
                'GNIP' : 'RASH',
                'CMUR' : 'RASH',
                'CPAI' : 'RASH',
                'TGON' : 'RASH',
                'EACE' : 'RASH',
                'CVEL' : 'RASH',
                'VBRA' : 'RASH',
                'PTET' : 'RASH',
                'TTHE' : 'RASH',
                'IMUL' : 'RASH',
                'PPER' : 'RASH',
                'SLEM' : 'RASH',
                'OTRI' : 'RASH',
                'SCOE' : 'RASH',
                'BNAT' : 'RASH',
                'PBRA' : 'RASH',
                'RFIL' : 'RASH',
                'MONO' : 'Excavata',
                'TVAG' : 'Excavata',
                'GINT' : 'Excavata',
                'SSAL' : 'Excavata',
                'LMAJ' : 'Excavata',
                'TBRU' : 'Excavata',
                'EGRA' : 'Excavata',
                'ADEA' : 'Excavata',
                'SCUL' : 'Excavata',
                'PHYT' : 'Excavata',
                'PERK' : 'Excavata',
                'BSAL' : 'Excavata',
                'NGRU' : 'Excavata',
                'GTHE' : 'Archaeplastida',
                'EHUX' : 'RASH',
                'CTOB' : 'RASH',
                'CPAR' : 'Archaeplastida',
                'GSUL' : 'Archaeplastida',
                'CMER' : 'Archaeplastida',
                'PPUR' : 'Archaeplastida',
                'CCRI' : 'Archaeplastida',
                'CVAR' : 'Archaeplastida',
                'CSUB' : 'Archaeplastida',
                'CREI' : 'Archaeplastida',
                'MNEG' : 'Archaeplastida',
                'VCAR' : 'Archaeplastida',
                'MSPE' : 'Archaeplastida',
                'OLUC' : 'Archaeplastida',
                'KFLA' : 'Archaeplastida',
                'MPOL' : 'Archaeplastida',
                'PPAT' : 'Archaeplastida',
                'SFAL' : 'Archaeplastida',
                'SMOE' : 'Archaeplastida',
                'PABI' : 'Archaeplastida',
                'ATRI' : 'Archaeplastida',
                'OSAT' : 'Archaeplastida',
                'NNUC' : 'Archaeplastida',
                'ATHA' : 'Archaeplastida',
                'ACOE' : 'Archaeplastida'}
supergroups2 = {sp : ('Opimoda' if group in ('Amoebozoa', 'Obazoa') else 'Diphoda') for sp, group in supergroups5.items()}