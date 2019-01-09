
# coding: utf-8

# In[1]:


# There isn't a Shannon Entropy script out there that disregards dashes that might throw off SE calculations. This script
# will take in a fasta file (each sequence must be restrained to only one line)


# In[15]:


import pandas as pd
import math


# In[10]:


def DashAlpha():
    Array1 = ["F", "L", "I", "M", "V","S", "P", "T", "A", "Y", "H", "Q", "N", "K", "D", "E", "C", "W", "R", "G","-"]
    return Array1


# In[11]:


def PureAlpha():
    Array1 = ["F", "L", "I", "M", "V","S", "P", "T", "A", "Y", "H", "Q", "N", "K", "D", "E", "C", "W", "R", "G"]
    return Array1


# In[129]:


def FirstGroupAlpha():
    Array1 = ['H','P','A','B','N']
    return Array1


# In[130]:


def SecondGroupAlpha():
    Array1 = ['H','P','A','B','N','S']
    return Array1


# In[13]:


def Sum(df, alphabet):
    TrueSum = 0
    for letter in alphabet:
        Locate = df.loc[df['Letter'] == letter]
        AddThis = Locate['Letter'].count()
        TrueSum += AddThis
    return TrueSum


# In[99]:


def ShannonEntropyCalculation(df,alphabet,totalletter):
    SEcalc = 0
    for letter in alphabet:
        LetterDF = df.loc[df['Letter'] == letter]
        LetterCount = LetterDF['Letter'].count()
        if (LetterCount != 0):
            Frequency = LetterCount/totalletter
            Calc = (-1 * Frequency) * (math.log2(Frequency))
            SEcalc += Calc
    return SEcalc
        


# In[81]:


def PosDFmaker(df,position):
    RowCount = 0
    PosDF = pd.DataFrame()
    Rows = df['Seq'].count()
    for x in range(0,Rows):
        Sequence = df.at[x,'Seq']
        PosDF.at[x,'Letter'] = Sequence[position]
    return PosDF


# In[38]:


def SeqDFmaker(file):
    Counter = 0
    Record = False
    SeqDF = pd.DataFrame()
    for line in open('FullSeqAlign.fa'):
        if(Record == True):
            SeqDF.at[Counter,'Title'] = Title
            SeqDF.at[Counter,'Seq'] = line
            Record = False
            Counter += 1
        if ('>' in line):
            Title = line[1:].strip()
            Record = True
    return SeqDF
    


# In[118]:


HGroup = ['G','A','L','V','I','M','P']
def ResiConvert(line,array,trueletter):
    for item in array:
        line = line.replace(item,trueletter)
    return line
        


# In[155]:


# First set refers to pure groupings based on side chains
def FirstSetGroupConversion(df):
    # Remember to use sequence!
    HGroup = ['G','A','L','V','I','M','P']
    PGroup = ['S','T','C','N','Q']
    AGroup = ['F','Y','W']
    BGroup = ['H','K','R']
    NGroup = ['D','E']
    OneGroup = ['1']
    TwoGroup = ['2']
    ThreeGroup = ['3']
    FourGroup = ['4']
    FiveGroup = ['5']
    RowCount = df['Seq'].count()
    for i in range(0,RowCount):
        TruthI = i
        OrigSeq = df.at[i,'Seq']
        Seq = ResiConvert(OrigSeq,HGroup,'1')
        Seq = ResiConvert(Seq,PGroup,'2')
        Seq = ResiConvert(Seq,AGroup,'3')
        Seq = ResiConvert(Seq,BGroup,'4')
        Seq = ResiConvert(Seq,NGroup,'5')
        Seq = ResiConvert(Seq,OneGroup,'H')
        Seq = ResiConvert(Seq,TwoGroup,'P')
        Seq = ResiConvert(Seq,ThreeGroup,'A')
        Seq = ResiConvert(Seq,FourGroup,'B')
        Seq = ResiConvert(Seq,FiveGroup,'N')
        df.at[TruthI,'Seq'] = Seq


# In[156]:


# Second set puts glycine and proline in their own spacces, specifically interested in 
# Amino acid structures. Will become handy if you have a lot of loops/coils
def SecondSetGroupConversion(df):
    # Remember to use sequence!
    HGroup = ['A','L','V','I','M']
    PGroup = ['S','T','C','N','Q']
    AGroup = ['F','Y','W']
    BGroup = ['H','K','R']
    NGroup = ['D','E']
    SGroup = ['G','P']
    OneGroup = ['1']
    TwoGroup = ['2']
    ThreeGroup = ['3']
    FourGroup = ['4']
    FiveGroup = ['5']
    SixGroup = ['6']
    RowCount = df['Seq'].count()
    for i in range(0,RowCount):
        TruthI = i
        OrigSeq = df.at[i,'Seq']
        Seq = ResiConvert(OrigSeq,HGroup,'1')
        Seq = ResiConvert(Seq,PGroup,'2')
        Seq = ResiConvert(Seq,AGroup,'3')
        Seq = ResiConvert(Seq,BGroup,'4')
        Seq = ResiConvert(Seq,NGroup,'5')
        Seq = ResiConvert(Seq,SGroup,'6')
        Seq = ResiConvert(Seq,OneGroup,'H')
        Seq = ResiConvert(Seq,TwoGroup,'P')
        Seq = ResiConvert(Seq,ThreeGroup,'A')
        Seq = ResiConvert(Seq,FourGroup,'B')
        Seq = ResiConvert(Seq,FiveGroup,'N')
        Seq = ResiConvert(Seq,SixGroup,'S')
        df.at[TruthI,'Seq'] = Seq
        return df
    


# In[161]:


AlphaAnswer = input("Hello! This is a Shannon Entropy script. Would you like to consider dashes in your script or not? (y/n)")
GroupAnswer = input("This script can also calculate SE based on amino acid groups. If you are not interested, enter n. If you are interested, the groupings are separated based on glycine/proline. Would you like these considered as hydrophobic -enter 1- or would you like them considered in their own group -enter 2-  ")
SequenceDF = SeqDFmaker('FullSeqAlign.fa')
if (AlphaAnswer == 'y'):
    UAlpha = DashAlpha()
if(AlphaAnswer == 'n'):
    UAlpha = PureAlpha()
if(AlphaAnswer.lower() != 'y' and AlphaAnswer.lower() != 'n'):
    print(AlphaAnswer.lower())
    print('Failure due to unreadable answer. Please restart script and try again. Please only use characters y or n for answers.')
if(GroupAnswer == '1'):
    FirstSetGroupConversion(SequenceDF)
    UAlpha = FirstGroupAlpha()
if(GroupAnswer == '2'):
    SecondSetGroupConversion(SequenceDF)
    UAlpha = SecondGroupAlpha()
# Code above can be implemeneted in script to ask for user input on whether or not to include dashes in calculations as '21st' aa

ResultSheet = pd.DataFrame()
SeqLength = len(SequenceDF.at[1,'Seq'].strip())
FirstSequence = SequenceDF.at[0,'Seq']
for i in range(0,SeqLength):
    TrueI = i
    PositionDF = PosDFmaker(SequenceDF,TrueI)
    SumValue = Sum(PositionDF,UAlpha)
    SEValue = ShannonEntropyCalculation(PositionDF,UAlpha,SumValue)
    ResultSheet.at[TrueI,'SeqPos'] = int(TrueI + 1)
    ResultSheet.at[TrueI,'FirstSequenceInput'] = FirstSequence[TrueI]
    ResultSheet.at[TrueI,'EntropyValue'] = SEValue

