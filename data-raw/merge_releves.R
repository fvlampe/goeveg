## Other function: Merge different releves
# merge_releves()

## Format: Vegetation table with samples in rows and species in columns
## How to define?
## Define matrix and define merging column (2 arguments)
## Based on definition in header data: one entry per sample new names for releves (e.g. water bodies, regions)
## How to aggregate? Methods:
## Presence/absence (simple/default);
## Transformation of "from"-scale to percentage values and calculation of mean (decide: treat zero values as NA or 0; default: 0) --> back-transformation in defined "to" scale
## "from" scale can be defined for all releves ("braun.blanquet", etc.) or as header column with definition (identical names) for each releve


## How to calculate?
# For samples in rows, the row sums must be calculated for common species;
## same algorithm as for merge_taxa using long table?
## Need test dataset!

#Bei merge releves:
# Releves müssen verschiedene Skalen haben können, basierend auf einer Scale Spalte (zusätzlich zur merge Spalte)
# (Cov2per)
# Bzw. Kann cov2per auf verschiedene teildatensätze mit verschiedenen Skalen angewendet werden !!
# (vermutlich einfacher und cov2per muss nicht verschiedene Teildatensätze können)
# merge_taxa wäre sicherlich gut über verschiedene Skalen; andererseits könnte

## VERGLEICH SYNTABLE, type = frequency or mean percentage ... --> im Prinzip gleiche Funktion
# passe Funktionen und Beschreibungen ggf. an syntable an

vegtable <- schedenveg
head <- schedenenv

# Testidee
# Fasse alle Aufnahmen in 2 Aufnahmen zusammen nach community
head$comm[head$comm == unique(head$comm)[1]]

# Make selection based on header column
unique(head$comm)[1]
unique(head$comm)[2]

sel <- schedenveg[head$comm == unique(head$comm)[2], ]

# Calculate new cover values for merged result
apply(sel > 0, 2, sum)  # Number of species in selected releves  (can be made as PA or proportion/frequency)
apply(sel, 2, sum)      # Sum of all coverages (should made relative, not very useful)
apply(sel, 2, mean)     # Mean of all coverages, could be nice

# Store merged result under name of group

stack(schedenveg)




