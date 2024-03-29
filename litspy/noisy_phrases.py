# Lists of noise to be used to eliminate noisy synonyms from synonym lists

# Some ontology nodes contain synonym lists in comment fields that can contain other information such as links,
# email addresses, citations, curation notes, definitions etc.
# Therefore, any collected terms that contain any of the following noise indicators should be removed

synonym_noise_indicators = [":", "@", " email ", "doi.org", "Wikipedia", "github", "TODO ", " et al", "th ed.", "[WP]",
                            "see also", "see article", "Editor node", "Editor note", "Taxon notes ",
                            "Consider merging", "mapping confirmed", "partof ", "Requires expert input",
                            "UMLS CUI", "synonyms", " doid ", "doid/", "Xref ",
                            "Definition based on", "characterized by", "symptoms ",
                            "believed to be derived from", "We place ",
                            "mice have ", "mouse has ", " to form ", "will be ceded", "use the term", "same name",
                            "presumed but not proven", "occurs in", "are different", "term renamed"]


# Terms that are not frequently-used synonyms for genes but that are frequently used in publications should be removed
# from gene synonym lists to increase precision (potentially at the cost of some recall). Examples include:
# three-letter amino acid codes (ASP),
# elements/compounds (CO2),
# measurement terms (mcL)
# month abbreviations
# biomedical words and abbreviations (RR, BP, CAS (in crispr-cas), traits, FDA),
# gene synonyms that are very common (e.g. xBP for various binding proteins results in many incorrect hits for all xBPs)
# statistical terms (cox regression model, r2 (r squared)) and other words/abbreviations common in publications (FIG)
# gene synonyms that are also words in various languages (AGE, ERA, ES, DANCE),


common_gene_noise = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                     "PRO", "PYL", "SEC", "SER", "THR", "TRP", "TYR", "VAL",
                     "CO2", "CO 2", "HCN", "MCL",
                     "JAN", "FEB", "MAR", "APR", "APRIL", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC",
                     "A11", "ABBA", "ABL", "ABP", "ABS", "ACCA", "ACK", "ACS", "ACT", "AD3", "AD4", "AD5", "ADP", "AFT",
                     "AGAS", "AGT", "AIK", "AIS", "AIT", "AKA", "ALBA", "ALC", "ALF", "ALP", "ALP", "ALS", "ALY",
                     "AMIGO", "AML 2", "AMP", "AMY", "ANA", "ANOVA", "APOLLO", "APPD", "APPL", "APX", "ARF", "ARIA",
                     "ARIA", "ARTEMIS", "ARX", "ASC", "ASPS", "AST", "ATAR", "AURA", "B10", "BACH", "BAFF", "BAL",
                     "BAM", "BAP", "BAP", "BEF", "BEN", "BENE", "BEY", "BIM", "BKS", "BLAST", "BLS", "BOCA", "BOD",
                     "BOM", "BOO", "BOP", "BOR", "BRL", "BTL", "BUN", "CAD", "CAD", "CAD", "CAL", "CALC", "CALIF",
                     "CALP", "CAM", "CAMP", "CANION", "CAP", "CAPB", "CAPER", "CAPON", "CAPRICE", "CARF", "CARF",
                     "CARK", "CAS", "CASPER", "CAV", "CCT", "CDF", "CDR", "CEE", "CERT", "CHA", "CHASM", "CHICA", "CHICO",
                     "CHIP", "CHIT", "CHN", "CHS", "CIA", "CIP", "CIR", "CIS", "CLAN", "CLAP", "CLASP", "CML", "COCO",
                     "COCOA", "CPD", "CPL", "CPU", "CRES", "CRIP", "CRL", "CSS", "CST", "CST", "CT2", "CTF", "CTR",
                     "CTS", "CYC", "D10", "DAG", "DALI", "DAMS", "DAN", "DANTE", "DAO", "DAP", "DAPPER", "DAT", "DBL",
                     "DBP", "DEF", "DEG", "DELE", "DELTA", "DENTS", "DHP", "DIA", "DICER", "DIF", "DIP", "DISP", "DIVA",
                     "DOM", "DOR", "DORA", "DOS", "DPS", "DRAGON", "DRAM", "DRT", "DTD", "DTP", "DUP", "EAD", "EDH",
                     "EEN", "EGO", "EKG", "ELKS", "EMP", "ENED", "ENGL", "ENIGMA", "ENL", "EOS", "EOS", "EPA", "EPI",
                     "EPIL", "ERIC", "ERIS", "ESP", "EST", "EXP", "FAB", "FAC", "FAD", "FAE", "FAG", "FAS", "FDA",
                     "FELS", "FIAT", "FIB", "FIP", "FIR", "FLATTOP", "FLT", "FOG", "FON", "FOP", "FPS", "FVL", "GADS",
                     "GAJ", "GALA", "GAT", "GATA", "GDS", "GEM", "GEN", "GENESIS", "GGR", "GIP", "GIT", "GLI", "GLOB",
                     "GOLIATH", "GOOFY", "GOR", "GOX", "GPCR", "GPD", "GPI", "GRAF", "GRAIL", "GRIPE", "GROG", "GRS", "GRX",
                     "GUP", "HAF", "HAI", "HAK", "HALP", "HAP", "HARE", "HARP", "HBP", "HCA", "HCC", "HDR", "HDRS",
                     "HED", "HEP", "HEP", "HERMES", "HERP", "HES", "HET", "HEX", "HIC", "HILI", "HIPPI", "HIR", "HOGA",
                     "HRG", "HRS", "HSM", "HYD", "HYL", "HYPE", "HYPERION", "ICEBERG", "IDOL", "IDP", "IFF", "IFI", "IMP",
                     "INF", "IPL", "IPS", "JAMA", "JAMB", "JUNO", "KAB", "KAF", "KAL", "KALI", "KAP", "KAT", "KEN",
                     "KEPI", "KET", "KGF", "KINO", "KIP", "KIST", "KLIP", "KOP", "KOR", "LAB", "LACS", "LAD", "LAH", "LAK",
                     "LAK", "LAN", "LAP", "LARGEN", "LAS", "LAT", "LAX", "LBP", "LCA", "LCA", "LECT", "LED", "LIB",
                     "LIM", "LIND", "LIPA", "LIR", "LIT", "LOR", "LPD", "LSK", "LUST", "LYRIC", "MACH", "MAD", "MAG",
                     "MAGMAS", "MAIL", "MAIR", "MAL", "MAL", "MANI", "MARC", "MARE", "MASA", "MAST", "MATER", "MCP",
                     "MCT", "MED", "MENT", "MEP", "MER", "MFR", "MGR", "MIB", "MIDAS", "MIM", "MIMA", "MINERVA",
                     "MINION", "MIR", "MIRK", "MIS", "MISE", "MMR", "MOCA", "MOLT", "MONA", "MONAD", "MOS", "MRS",
                     "MSF", "MSS", "MTC", "MTD", "MTS", "MUD", "MUSTANG", "MUT", "MYG", "MYM", "NAF", "NAF", "NAG",
                     "NAK", "NAM", "NAN", "NAP", "NAPA", "NAR", "NARR", "NAT", "NBS", "NDF", "NEMO", "NEP", "NESH",
                     "NIP", "NIPA", "NIS", "NIX", "NKR", "NOBODY", "NOS", "NOXA", "NPI", "NUANCE", "O11", "OASIS",
                     "OBOE", "OPS", "OPT", "ORF", "ORF", "OSSA", "OVAL", "PACT", "PAD", "PAL", "PAL", "PAM", "PAP",
                     "PAPA", "PAPAS", "PAR", "PARC", "PARI", "PARIS", "PARS", "PATE", "PAUL", "PBS", "PCP 2", "PCR",
                     "PED", "PEGASUS", "PENUMBRA", "PEP", "PEP", "PEPS", "PERF", "PES", "PESKY", "PICH", "PICOT",
                     "PIKA", "PILAR", "PIPPIN", "PIS", "PIST", "PKG", "PLAP", "PMK", "POLK", "POTE", "PPD", "PPH",
                     "PPH", "PPT", "PRAT", "PREP", "PRIMA", "PRISM", "PRN", "PRP", "PSF", "PSST", "PST", "PST", "PTA",
                     "PTC", "PTG", "PTP", "PTP", "PURL", "RAC1", "RAD", "RAGA", "RAH", "RAMP", "RASI", "RAX", "REA",
                     "REC", "REGR", "REN", "RHA", "RHOS", "RHS", "RIFF", "RISC", "RISP", "RIT", "RITA", "RNS", "ROG",
                     "ROM", "ROS", "ROS", "ROX", "RSS", "SAA", "SAB", "SAC", "SAG", "SAGE", "SAHH", "SALSA", "SAN",
                     "SANCHO", "SANS", "SAP", "SAP", "SAPS", "SARI", "SCAD", "SCAP", "SCF", "SCH", "SCOP", "SCOT", "SDS",
                     "SECT", "SELS", "SELT", "SEME", "SEP", "SEP", "SERA", "SERS", "SGD", "SHANK", "SHAPY", "SHP",
                     "SIKE", "SIL", "SIMP", "SISE", "SIVA", "SLA", "SLAT", "SLD", "SLICK", "SLOB", "SLT", "SLY", "SMIT",
                     "SNARK", "SOLO", "SONE", "SONE", "SOUL", "SPAK", "SPARTAN", "SPL", "SPP", "SPR", "SPS", "SPS",
                     "SPT", "SRA", "STA", "STD", "STG", "STP", "STR", "STRAD", "STS", "SUP", "SWA", "SYL", "SYM", "SYN",
                     "TAJ", "TALI", "TANGO", "TAPA", "TARA", "TAU", "TC1", "TCB", "TCI", "TEAP", "TECH", "TECK", "TED",
                     "TEL", "TELE", "TEM", "TER", "TERA", "TERP", "TES", "TGT", "TIAR", "TIC", "TKT", "TLN", "TMC", "TMS",
                     "TNT", "TOB", "TOB", "TOM", "TOR", "TRAD", "TRAG", "TRF", "TRID", "TRP", "TRT", "TSK", "TSP",
                     "TUBA", "TULA", "TYP", "UFO", "UGT", "UNRIP", "URB", "UTI", "VAN", "VASA", "VEL", "VIII", "VIN",
                     "VIP", "VISTA", "WABS", "WARP", "WBS", "WICH", "WID", "YAP", "YETI", "YRS", "ZAC", "ZAG", "ZAK",
                     "ZAP",
                     "AMINO ACID TRANSPORTER", "BINDING PEPTIDE", "HYDROLASE", "PORIN", "PROTEIN C",
                     "RNA PROCESSING FACTOR",
                     "ELK", "ERB", "ERK", "GBP", "MST", "MPP", "P24", "P25", "P35", "P36", "P38", "P57", "P75", "P100",
                     "P200", "RAB", "RBP",
                     "AIM 1", "AIM 2", "AIM", "COX", "D10", "EPO", "FIG", "LCA", "PCA", "REF", "TOP 2",
                     "ACT", "AFAR", "AFRO", "AGE", "AGO", "AID", "AIR", "ALIEN", "ALL", "APE", "APP", "APPS", "APT",
                     "ARC", "ARC", "ARCH", "ARK", "ARM", "ARMER", "ARMS", "ART", "ARTS", "ASAP", "ASK", "ATOPY", "AURA",
                     "BANK", "BAR", "BARS", "BART", "BASE", "BASH", "BAT", "BEST", "BIKE", "BIT", "BITE",
                     "BLAME", "BOG", "BOMB", "BOR", "BOULE", "BRAG", "BRAVO", "BRIGHT", "BRUCE", "CAGE", "CAIN", "CALM",
                     "CAMEL", "CAN", "CAP", "CAPS", "CAR", "CARDINAL", "CARMEN", "CARP", "CARP", "CART", "CASH", "CAST",
                     "CATS", "CAVA", "CHAMP", "CHIMP", "CHOP", "CIG", "CINEMA", "CLAMP", "CLINT", "CLIP", "COASTER",
                     "COD", "COP", "COT", "CRAM", "CRAMP", "CREPT", "CREST", "CROP", "CUT", "DAMAGE", "DANCE", "DANGER",
                     "DEAR", "DEEPEST", "DEFT", "DES", "DIETER", "DINE", "DING", "DREAM", "EAR", "EARS", "END", "ENRAGE",
                     "ERA", "FACT", "FAD", "FAME", "FAN", "FAST", "FAT", "FATE", "FATS", "FELL", "FETA", "FIND", "FISH",
                     "FIX", "FIX", "FLAME", "FLAP", "FLASH", "FLIP", "FOE", "FOG", "FOR", "FRA", "FRITZ", "GAP", "GAS",
                     "GASP", "GET", "GIF", "GILT", "GOA", "GOAT", "GRAB", "GREAT", "GRIT", "GULP", "HAD", "HANK",
                     "HARP", "HASNT", "HEED", "HELIOS", "HIP", "HITS", "HOP", "HUB", "HUG", "ICE", "INCA", "INCL",
                     "IOTA", "JAB", "KID", "KILLER", "LAG", "LAMP", "LAP", "LARD", "LARK", "LES", "LETS",
                     "LIAR", "LIFEGUARD", "LIGHT", "LIME", "LIP", "LIT", "LOBE", "LORD", "LUST", "MAI", "MARK", "MART",
                     "MASK", "MASK", "MASS", "MAT", "MATT", "MATTER", "MEMO", "MEN", "MES", "MICE", "MINK", "MINK",
                     "MINOR", "MINT", "MIST", "MOB", "MOP", "MORT", "NAIL", "NEST", "NET", "NET", "NETS", "NEU", "NOPE",
                     "NOT", "NUDE", "NUT", "ODD", "ORCA", "OUT", "PACER", "PALLID", "PANDA", "PANDER", "PARTICLE", "PEAS", "PEN",
                     "PERK", "PILOT", "PILOT", "PIN", "PIN", "PINCH", "PINS", "PINT", "PLEIAD", "PLUTO", "POEM", "POSH",
                     "POSHER", "PREY", "PUMA", "PUNISHER", "RACE", "RAGS", "RAIN", "RAM", "RAMP", "RANK", "RAY", "RED",
                     "RHINO", "RHO", "RICK", "RIG", "RIM", "RIM", "RIP", "RIP", "ROD", "SANG", "SCAR", "SCRAPS",
                     "SECRET", "SEX", "SHARP", "SHIP", "SHOT", "SIMPLE", "SIN", "SIP", "SIT", "SKIP", "SKY", "SLACK",
                     "SLAP", "SMILE", "SNAIL", "SNIP", "SPAR", "SPASM", "SPICE", "SPIN", "SPRIGHTLY", "SPRING", "STAR",
                     "STARING", "STARS", "STELLAR", "STEP", "STING", "STOP", "STRAP", "STUD", "SWAN", "SWAP", "TACTILE",
                     "TAG", "TAP", "TAP", "TASK", "TAPS", "TAUT", "TEMP", "THANK", "THE", "THETA", "TIED", "TIM", "TIM",
                     "TIM", "TIP", "TOP", "TRADE", "TRAIL", "TRAIL", "TRAITS", "TRAM", "TRAMP", "TRANCE", "TRAP",
                     "TRAP", "TRIM", "TRIP", "TROY", "TUBE", "TUG", "TUG", "TUNA", "TWEAK", "TWINKLE", "TYPE", "VISA",
                     "WAR", "WARS", "WARTS", "WAS", "WASP", "WASPS", "WAVE", "WAVE", "WHIP", "WHISTLE", "WIRE", "WISH",
                     "WISP", "YES", "ZETA", "ZIP"
                     ]

stop_words = ["a", "able", "about", "across", "after", "all", "almost", "also", "am", "among", "an", "and", "any",
              "are", "as", "at", "be", "because", "been", "but", "by", "can", "cannot", "could", "dear", "did", "do",
              "does", "either", "else", "ever", "every", "for", "from", "get", "got", "had", "has", "have", "he",
              "her", "hers", "him", "his", "how", "however", "i", "if", "in", "into", "is", "it", "its", "just",
              "least", "let", "like", "likely", "may", "me", "might", "most", "must", "my", "neither", "no", "nor",
              "not", "of", "off", "often", "on", "only", "or", "other", "our", "own", "rather", "said", "say", "says",
              "she", "should", "since", "so", "some", "than", "that", "the", "their", "them", "then", "there", "these",
              "they", "this", "tis", "to", "too", "twas", "us", "wants", "was", "we", "were", "what", "when",
              "where", "which", "while", "who", "whom", "why", "will", "with", "would", "yet", "you", "your"]

not_top_ten = ["de", "la", "el", "en", "lo", "del", "que", "et", "le", "un", "du", "des",
               "medical", "clinical", "health", "hospital", "participant", "participants", "patient", "patients",
               "scheme", "schemes", "program", "programs", "programme", "programmes", "center", "centers", "centre",
               "centres",
               "abstract", "introduction", "background", "method", "methods", "result", "results", "conclusion",
               "conclusions", "findings", "outcome", "evidence", "aim", "aims", "doi",
               "gene", "genes", "genome", "proteins", "protein", "cell", "cells", "disease", "diseases", "syndrome",
               "syndromes", "symptom", "symptoms",
               "test", "tests", "study", "studies", "work", "experiment", "experiments", "technique", "techniques",
               "technical", "analysis", "analyses", "analysing", "criteria",
               "day", "days", "time", "hour", "hours", "recent", "during", "effect", "effects", "used", "using",
               "compare", "compared", "comparison", "confirm", "confirmed", "confirmation",
               "association", "associated", "available", "availability",
               "previously", "report", "reports", "reported", "unreported",
               "many", "more", "less", "fewer", "number", "numbers", "large", "larger",
               "level", "levels", "total", "different", "further", "such",
               "well", "including", "being", "within", "anti", "data", "show", "shown",
               "case", "cases", "control", "controls",
               "one", "two", "three", "four", "five"
               ]
