library(randomForest)

# load the model
load('/Users/chengwang/Dropbox/ws/dvdep/deltavina/deltavina/models/rffit.rda')

# input and output file name
infn = 'input.csv'
outfn = 'output.csv'

print(paste("Read input: ", infn))
# read in input as dataframe df
df = read.table(infn, header=T, stringsAsFactors = F, sep=',')

print(df)

# get features from df
feats = df[3:22]

# predict the binding affinity
pred = round(predict(rffit, newdata = feats) + df$vina,2)

# write output
output = data.frame(pdb = df$pdb, pred = pred)

print(paste("Write input: ", outfn))
write.table(output, outfn, sep=',', row.names = F, quote = F)

print("Done")
    
