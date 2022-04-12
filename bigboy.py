import csv
from csv import reader
import copy
import sys
import re
import pandas as pd
from pandas import *
from Bio.Align import PairwiseAligner
from colour import Color
from pprint import pprint
import datetime
#opens file required for reading
def openFile(fileName):
    data = open(fileName+'.csv')
    rows = []
    for row in reader(data):
        rows.append(row)
    data.close()
    return rows

def closeAndSave(fileName, dataBase):
    MONTH,DATE,HOUR,MINUTE = str(datetime.date.today().month),str(datetime.date.today().day),str(datetime.datetime.now().hour),str(datetime.datetime.now().minute)

    timeStamp = "-FIX_"+DATE+'of'+MONTH+'at'+HOUR+MINUTE+'.csv'
    with open(fileName+timeStamp, 'w') as f:
        write = csv.writer(f)
        for row in range(0,len(dataBase)):
            write.writerow(dataBase[row])
        print('File Saved As:'+fileName+timeStamp)

def squiggly(value, custom = None):#display purposes, send a string to be surrounded. custom is for custom squiqqles 
    squiggle = '~@~@~@~@~@~@~@~@~@~'

    if(custom is None):
        print(squiggle+'|'+value+'|'+squiggle+'\n')
    else:
        print(squiggle+'|'+custom+'|'+squiggle+'\n')

def displayTable(rows, fileName):#outputs the column headers prettily
    top10 = []
    squiggly(fileName)
    for row in range(0,10):
        top10.append(rows[row])
    print(DataFrame(top10))

def deleteColumn(rows):
    columnIndex = int(input("which column do you want to be removed? enter index: "))
    removed = str(rows[0][columnIndex])
    for row in rows:
        row.pop(columnIndex)
    return rows

def insertColumn(rows, columnIndex, columnName):
    rows[0].insert(columnIndex,columnName)
    for row in range(1,len(rows)):
        rows[row].insert(columnIndex,'')
    return rows

def checkMembership(colourString,colourRepo,flag):
    score,tempScore,aligner,bestColour = 0,0,PairwiseAligner(),''

    if(colourString == flag):
        return flag
    
    words = ((((str(colourString).strip()).lower()).replace('/',' ')).replace('melange','')).split(' ')#remove melange, its a fabric type not colour related
    for element in range(0,len(words)):
        word = str(words[element]).strip()
        if(word.lower() == 'melange'):
            continue
        if word in colourRepo:
            return word
        else:
            for colour in colourRepo:
                if(len(word)==0):
                    continue
                else:
                    tempScore = aligner.score(word,colour)
                    if score < tempScore:
                        score = tempScore
                        bestColour = colour
    return bestColour

def colourArray(rows,colourSource,colourDestination,flag):
    line,preservedColours,skip1st = '',[],0

    with open('assets/colours.txt') as f:
        colours = f.readlines()

    for colour in range(0,len(colours)):
        line = colours[colour]
        line = (line.strip()).lower()
        tempPair = line.split(',')
        preservedColours.append(tempPair)
        line = tempPair[0]
        colours[colour] = line

    for row in rows:
        if(skip1st == 0):
            skip1st = 1
            continue
        colourLeftPair = checkMembership(str(row[colourSource]),colours,flag)
        for pair in preservedColours:
            if(colourLeftPair.lower()==(pair[0]).lower()):
                row[colourDestination] = pair[1]
                continue
    return rows

def generateParents(tableC,rangeIdentifier,columnsToKeep):
    parentColumn,Insert1st,ranInt,tableP,currentSKUrange, = rangeIdentifier+1,0,0,[],""
    tableC = insertColumn(tableC, parentColumn, 'Parent SKU')#inserts parent sku column
    tableP.append(tableC[0])#add headers to new list that we will build on

    for row in range(1,len(tableC)):
        ranInt+=1
        if(currentSKUrange!=tableC[row][rangeIdentifier]):#new range if true
            ranInt,currentSKUrange,tempclone = 0,tableC[row][rangeIdentifier],copy.deepcopy(tableC[row])

            temp = [''] * len(tableC[row])
            for column in range(0,len(columnsToKeep)):#gets columns to carry to parent
                index = int(columnsToKeep[column])
                temp[index] = tableC[row][index]
            temp[rangeIdentifier] = currentSKUrange
            tableP.append(temp)#parent

            tempclone[rangeIdentifier] = str(currentSKUrange)+str(ranInt)
            tempclone[parentColumn] = currentSKUrange
            tableP.append(tempclone)#firstchild
        else:#same range
            tableC[row][rangeIdentifier] = str(currentSKUrange)+str(ranInt)
            tableC[row][parentColumn] = currentSKUrange
            tableP.append(tableC[row])#child
    return tableP

def collectValuesForParent(rows,targetColumnIndex,parentColumnIndex,delim):
    values,currentParentIndex,currentSKUrange = [],0,""

    for row in range(1,len(rows)):
        if(rows[row][targetColumnIndex].strip() == ""):
            rows[currentParentIndex][targetColumnIndex] = delim.join(list(set(values)))
            values,currentParentIndex = [],int(row)
        else:#same range    
            values.append(str(rows[row][targetColumnIndex]).strip())
    return rows

def customFunction(rows,columnIndex):
    functions = ['multiply by','divide by','to x add','to x minus']
    for function in range(0,len(functions)):
        print(function,": ",functions[function])
    choice = input("What would you like to do to "+rows[0][columnIndex]+"? choose by index, then value with n to seperate e.g. 2n44.2 = add 44.2: ")
    choice = choice.split('n')
    functionChoice = int(choice[0])
    operand = float(choice[1])
    if(functionChoice==0):
        for row in range(1,len(rows)):
            if(str(rows[row][columnIndex])!='' or str(rows[row][columnIndex])!=' '):
                rows[row][columnIndex] = float(rows[row][columnIndex])*operand
    elif(functionChoice==1):
        for row in range(1,len(rows)):
            if(str(rows[row][columnIndex])!='' or str(rows[row][columnIndex])!=' '):
                rows[row][columnIndex] = float(rows[row][columnIndex])/operand
    elif(functionChoice==2):
        for row in range(1,len(rows)):
            if(str(rows[row][columnIndex])!='' or str(rows[row][columnIndex])!=' '):
                rows[row][columnIndex] = float(rows[row][columnIndex])+operand
    elif(functionChoice==3):
        for row in range(1,len(rows)):
            if(str(rows[row][columnIndex])!='' or str(rows[row][columnIndex])!=' '):
                rows[row][columnIndex] = float(rows[row][columnIndex])-operand
    return rows

def mergeColumns(rows, column1Index, column2Index, targetColumn, newColumnFlag, delim):
    if not newColumnFlag:
        rows = insertColumn(rows, targetColumn, str(rows[0][column1Index])+' + '+str(rows[0][column2Index]))
    for row in rows:
        row[targetColumn] = str(row[column1Index])+delim+str(row[column2Index])
    return rows

functions,delim = ['Generate Parents','Generate Generic Colours','Collate Values Into Parents','Delete Column','Custom Transformation','Merge Columns','Swap Columns - Placeholder','Close and Save'],'|'
filename = input("path(if not in current dir) and name of csv file excluding type e.g. 'myexports/export1': ")
rows = openFile(filename)
print("@~@~ warning - all indexes are counted from 0. e.g. the first column is 0 second is 1 ~@~@")
print("@~@~ warning - Parent generation must be performed first ~@~@")
while True:
    print("Hi, this program generates parent/child relationships among other functions: ")
    for function in range(0,len(functions)):
        print(function,": ",functions[function])
    functionChoices = input("list the functions you want to be performed e.g. '0,2,3': ").strip()
    functionChoices = functionChoices.split(',')
    functionChoices = [int(i) for i in functionChoices]

    if 0 in functionChoices:
        squiggly('Generate Parents')
        rangeIdentifier = int(input('Which column can be used to identify variation ranges? e.g. item group, barcode. enter as index: '))
        tempColumnsToKeep = input('which columns would you like to be preserved in the new parent rows? enter as list of indexes e.g. 5,6,7: ').strip()
        tempColumnsToKeep = tempColumnsToKeep.split(',')
        columnsToKeep = [int(i) for i in tempColumnsToKeep]
        rows = generateParents(rows,rangeIdentifier,columnsToKeep)
        displayTable(rows, filename)
        squiggly('')
        squiggly('complete')
        squiggly('')
    if 1 in functionChoices:
        squiggly('Generate Generic Colours')
        colourSource = int(input('Which column has the source colour? enter as index: '))
        flag = input('enter values you want to be skipped. e.g. \'no colour\': ')
        choice = input('do you want to insert a column to the right? if not they will overwrite the current colour column. y/n: ')
        if(choice.lower() == 'y'):
            rows = insertColumn(rows, colourSource+1, 'Generic Colour')
            rows = colourArray(rows,colourSource,colourSource+1,flag)
        elif(choice.lower() == 'n'):
            rows = colourArray(rows,colourSource,colourSource,flag)
        displayTable(rows, filename)
        squiggly('')
        squiggly('complete')
        squiggly('')
    if 2 in functionChoices:
        squiggly('Collate Values Into Parents')
        parentColumnIndex = int(input('Which column can be used to identify variation ranges? e.g. item group, barcode. enter as index: '))
        targetColumnIndex = int(input('Which columns values need to be collated into their parent row? enter as index: '))
        tempDelim = str(input('what delimiter would you like to be used for collated values? | is default for empty response: '))
        if(tempDelim!=''):
            rows = collectValuesForParent(rows,targetColumnIndex,parentColumnIndex,tempDelim)
        else:
            rows = collectValuesForParent(rows,targetColumnIndex,parentColumnIndex,delim)
        displayTable(rows, filename)
        squiggly('')
        squiggly('complete')
        squiggly('')
    if 3 in functionChoices:
        squiggly('Delete Column')
        rows = deleteColumn(rows)
        displayTable(rows, filename)
        squiggly('')
        squiggly('complete')
        squiggly('')
    if 4 in functionChoices:
        squiggly('Custom Transformation')
        columnIndex = int(input("Which column do you want to perform the function on? enter as index: "))
        rows = customFunction(rows,columnIndex)
        displayTable(rows, filename)
        squiggly('')
        squiggly('complete')
        squiggly('')
    if 5 in functionChoices:
        squiggly('Merge Columns')
        tempDelim = str(input('what delimiter would you like to be used for merged values? | is default for empty response: '))
        column1Index = int(input('column 1, enter as index: '))
        column2Index = int(input('column 2, enter as index: '))
        targetColumn = int(input('destination column, enter as index: '))
        tempFlag = input("Does the destination column already exist? y/n: ")
        if(tempFlag.lower() == 'y'):
            newColumnFlag = True
        elif(tempFlag.lower() == 'n'):
            newColumnFlag = False
        if(tempDelim!=''):
            rows = mergeColumns(rows, column1Index, column2Index, targetColumn, newColumnFlag, tempDelim)
        else:
            rows = mergeColumns(rows, column1Index, column2Index, targetColumn, newColumnFlag, delim)
        displayTable(rows, filename)
        squiggly('')
        squiggly('complete')
        squiggly('')
    if 6 in functionChoices:
        print('placeHolder')
    if 7 in functionChoices:
        closeAndSave(filename,rows)
        break
sys.exit()

#needs index of some value(column attribute) thats common between would be children to identify the range. e.g. product group, sku, barcode. tablep for new parents table, tablec is old table
#rangeIdentifier will be for skus, rangeIdentifier + 1 is where the parent SKU will be. user needs to make blank column to the right of their identifier


# filename = input("filename: ")
# cSource = input("colour source column: ")
# cDestination = input("new colour destination column: ")


# filename = "test files/baseCOT"
# cSource = 7
# cDestination = 16
# flag = "no colour"#for cells to be ignored
# rangeIdentifier = 0
# parentColumnIndex = rangeIdentifier+1#for collating values to parents, used to track parent range
# columnsToKeep = [4,5,6]


# def displayTable(rawList, fileName):#outputs the column headers prettily
#     currentIndex,whiteSpace = 0,'.'

#     squiggly(fileName)

#     print(end = '|')
#     for columnHeader in rawList[0]:#outputs indexing
#         length = round((len(columnHeader)-1)/2)
#         print(length*whiteSpace+str(currentIndex)+length*whiteSpace, end = '|')
#         currentIndex+=1
#     print('\n')

#     print(end = '|')#outputs values
#     for columnHeader in rawList[0]:
#         print(columnHeader, end = '|')
#     print('\n')


#this will only return the closest match to a knowable colour, to do: pair knowable colours to their generic countrerpart.
#however, hpoefulyl most of the time itll grab the generic one
#baseColours = ['red','green','blue','orange','yellow','black','purple','white']
# tempScore = aligner.score(stringInValue,colourCompare)
#     aligner,score,bestMatch,regex = PairwiseAligner(),0,'no match',re.compile('[^a-zA-Z ]')#colourRepo read as txt, each line is sky,blue. 0 = weird name 1 = simple name
#     baseColours = ['white','yellow','red','green','blue','orange','black','purple']
#     value = copy.deepcopy(rawValue).lower()
#     value = regex.sub('',value)


# fileName = 'test files/baseCot'
# closeAndSave(fileName,simpleColourController(openFile(fileName),8,9))


# with open('assets/colours.txt') as f:
#     colours = f.readlines()
# count = 0
# for line in colours:
#     tempPair = line.split(',')
#     line = str(tempPair[0]).lower()
#     print(count,':',line)
#     count+=1
# pprint(pg.color.THECOLORS)


# def colourArray(w3ColourReso,csColourReso,ogColourReso,rows,colourSource,colourDestination):
#     tempCompareSet,w3Best,csBest,ogBest,compareSet,w3csCheck,w3ogCheck,csogCheck,aligner = [],'','','',[],0,0,0,PairwiseAligner()
#     bestRows = ['']*(rows)
#     for row in range(0,len(rows)):
#         w3Best = checkMembership(copy.deepcopy(rows[row][colourSource]),w3ColourReso)
#         csBest = checkMembership(copy.deepcopy(rows[row][colourSource]),csColourReso)
#         ogBest = checkMembership(copy.deepcopy(rows[row][colourSource]),ogColourReso)
#         compareSet.append(w3Best)
#         compareSet.append(csBest)
#         compareSet.append(ogBest)
#         tempCompareSet = copy.deepcopy(compareSet)
#         if(len(set(compareSet))==1):
#             bestRows[row] = w3Best#append any of them, theyre the same value
#         else:
#             w3csCheck = aligner.score(w3Best,csBest)
#             w3ogCheck = aligner.score(w3Best,ogBest)
#             csogCheck = aligner.score(csBest,ogBest)
#             if(w3csCheck > w3ogCheck and w3csCheck > csogCheck):
#                 bestRows[row] = w3Best
#             elif(w3ogCheck > w3csCheck and w3ogCheck > csogCheck):
#                 bestRows[row] = ogBest
#             elif(csogCheck > w3ogCheck and csogCheck > w3csCheck):
#                 bestRows[row] = csBest


# def simpleColour(rawValue, buried, colourRepo):#identifies the base colour given a name or string, buried flag to indicate which
#     aligner,score,bestMatch,regex = PairwiseAligner(),0,'no match',re.compile('[^a-zA-Z ]')#colourRepo read as txt, each line is sky,blue. 0 = weird name 1 = simple name
#     baseColours = ['white','yellow','red','green','blue','orange','black','purple']
#     value = copy.deepcopy(rawValue).lower()
#     value = regex.sub('',value)
#     winner = ''
#     if(buried):
#         value = value.split(' ')
#     for colour in colourRepo:
#         tempColourValuePair = colour.split(',')
#         colourCompare = tempColourValuePair[0].lower()
#         if(buried):
#             for stringInValue in value:#splits string into words and checks each one
#                 if(stringInValue!='' and stringInValue !=' '):
#                     stringInValue = stringInValue.strip()
#                     tempScore = aligner.score(stringInValue,colourCompare)
#                     if(tempScore > score):
#                         bestMatch = tempColourValuePair[1] 
#                         score = tempScore
#                         winner = colourCompare
#                     elif(stringInValue in baseColours):
#                         return stringInValue,0,'x'
#         else:
#             if(value!='' and value !=' '):
#                 tempScore = aligner.score(value,colourCompare)
#                 if(tempScore > score):
#                     bestMatch = tempColourValuePair[1]
#                     score = tempScore
#                     winner = colourCompare
#                 elif(value in baseColours):
#                     return value,0,'x'
#     return bestMatch,score,winner


# def simpleColourController(rawList, colourIndex, simpleColourIndex):#feeds the row colour values into simplecolour
#     with open('assets/colours.txt','r') as file:
#         colours = file.readlines()
#     for row in range(1,len(rawList)):
#         rawList[row][9],rawList[row][10],rawList[row][11]  = simpleColour(rawList[row][colourIndex],False,colours)
#     return rawList


