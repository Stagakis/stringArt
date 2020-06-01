import numpy as np
import cv2
import math
from math import sin,cos,floor, sqrt
import random
import time

target_image_name = "target_final_new.png"
target_preprocessed = True

image_size = 640 #Assume square image
offset_x = image_size/2
offset_y = image_size/2

number_of_genes = 600
intial_gene_length = 400
survival_threshold = 0.90
mutate_chance = 0.05
grow_chance = 0.15

number_of_pins = 360
pin_size = 2
line_thickness = 1
line_color = 0

class Gene:
    valid_strings = np.zeros((number_of_pins,number_of_pins), np.uint8)

    def __init__(self, gene_length = -1):
        if(gene_length == -1):
            return
        self.pins = [random.randint(0, number_of_pins - 1) for i in range(0, gene_length)]
        self.cleanDuplicates()

    def cleanDuplicates(self):
        self.pins = [self.pins[i] for i in range(0, len(self.pins) - 1) if Gene.valid_strings[self.pins[i],self.pins[i+1]] ]
        self.pins.append(self.pins[-1])

    def mutate(self):
        index = random.randint(0, len(self.pins) - 1)
        while(True):
            new_pin = random.randint(0, number_of_pins - 1)
            if(Gene.valid_strings[new_pin,self.pins[index]] and Gene.valid_strings[new_pin,self.pins[max(0,index-1)]] and Gene.valid_strings[new_pin,self.pins[min(len(self.pins) - 1, index+1)]]):
                self.pins[index] = new_pin
                break
            else:
                continue

    def grow(self):
        while(True):
            new_pin = random.randint(0, number_of_pins - 1)
            if(Gene.valid_strings[new_pin,self.pins[-1]]):
                self.pins.append(new_pin)
                break
            else:
                continue

    def getGeneImage(self, pins_global):
        gene_image = np.full((image_size,image_size), 0, np.uint8)
        pin_indeces = self.pins
        [cv2.line(gene_image, (floor(offset_x + pins_global[pin_indeces[i]][0]), floor(offset_y +pins_global[pin_indeces[i]][1])), (floor(offset_x + pins_global[pin_indeces[i+1]][0]), floor(offset_y + pins_global[pin_indeces[i+1]][1])), 255-line_color, line_thickness) for i in range(0, len(pin_indeces) - 1)]
        return gene_image

def compareImages(gene_image, target):
    new_target = np.add(gene_image,target,dtype=np.int32)
    dist = np.count_nonzero(new_target==0)

    #a = np.array(new_target, dtype=np.uint8)
    #cv2.imshow("NEWTARGET", a)
    #cv2.waitKey(0)

    return np.sum(dist)

def evaluateGene(gene, target, pins):
    gene_image = gene.getGeneImage(pins)

    #gene_image = np.full((image_size,image_size), 0, np.uint8)
    #pin_indeces = gene.pins
    #[cv2.line(gene_image, (floor(offset_x + pins[pin_indeces[i]][0]), floor(offset_y + pins[pin_indeces[i]][1])), (floor(offset_x + pins[pin_indeces[i+1]][0]), floor(offset_y + pins[pin_indeces[i+1]][1])), 255-line_color, line_thickness) for i in range(0, len(pin_indeces) - 1)]

    length_factor = 1
    if(len(gene.pins) > intial_gene_length):
        length_factor = length_factor + 0.02*(len(gene.pins) - intial_gene_length)

    return compareImages(gene_image, target)*length_factor
    
def generatePins(N):
    radius = image_size/2
    step = (2*math.pi)/N
    pins = [[radius*cos(step*i), radius*sin(step*i)] for i in range(0, N)]
    return pins

def drawLine(canvas, pin1, pin2):
    point1 = (floor(offset_x + pin1[0]), floor(offset_y + pin1[1]))
    point2 = (floor(offset_x + pin2[0]), floor(offset_y + pin2[1]))
    cv2.line(canvas, point1, point2, line_color, line_thickness)

def drawGene(gene, pins): #OBSOLETE, TESTING PURPOSES ONLY
    gene_image = np.full((image_size,image_size), 255, np.uint8)
    pin_indeces = gene.pins
    [cv2.line(gene_image, (floor(offset_x + pins[pin_indeces[i]][0]), floor(offset_y + pins[pin_indeces[i]][1])), (floor(offset_x + pins[pin_indeces[i+1]][0]), floor(offset_y + pins[pin_indeces[i+1]][1])), line_color, line_thickness) for i in range(0, len(pin_indeces) - 1)]
    cv2.imshow("GEN", gene_image)

def saveGene(gene, pins, name, generation, score): #OBSOLETE, TESTING PURPOSES ONLY
    gene_image = np.full((image_size,image_size), 255, np.uint8)
    pin_indeces = gene.pins
    [cv2.line(gene_image, (floor(offset_x + pins[pin_indeces[i]][0]), floor(offset_y + pins[pin_indeces[i]][1])), (floor(offset_x + pins[pin_indeces[i+1]][0]), floor(offset_y + pins[pin_indeces[i+1]][1])), line_color, line_thickness) for i in range(0, len(pin_indeces) - 1)]
    cv2.imwrite("genes/GEN" + "_" + str(name) +"_" + str(generation) +"_"+ str(score)+ "_" + str(len(gene.pins)) +"_" + ".png", gene_image)

def drawPins(canvas, pins):
    for i in range(0, len(pins)):
        x_pos = floor(offset_x + pins[i][0])
        y_pos = floor(offset_y + pins[i][1])
        #print(str(x_pos) + "  " + str(y_pos))
        canvas[y_pos-pin_size:y_pos+pin_size, x_pos-pin_size:x_pos+pin_size] = 0

def createNewGene(gene1, gene2):
    
    new_gene = Gene()
    new_gene2 = Gene()
    #return new_gene, new_gene2

    length1 = len(gene1.pins)
    length2 = len(gene2.pins)
    min_length = min(length1, length2)

    point1 = random.randint(0, min_length - 10)
    point2 = random.randint(point1 + 1, min_length)

    new_gene.pins = gene1.pins.copy()
    new_gene2.pins = gene2.pins.copy()

    new_gene.pins[point1:point2]  = gene2.pins[point1:point2]
    new_gene2.pins[point1:point2] = gene1.pins[point1:point2]

    new_gene.cleanDuplicates()
    new_gene2.cleanDuplicates()

    new_gene.mutate()
    new_gene2.mutate()

    if(gene1.pins == new_gene.pins or gene1.pins == new_gene2.pins or gene2.pins == new_gene.pins or gene2.pins == new_gene2.pins):
        print("Point1/Point2 " + str(point1) +"/"+str(point2) )
        print("Old1 " + str(gene1.pins))
        print("Old2 " + str(gene2.pins))
        
        print("New1 " + str(new_gene.pins))
        print("New2 " + str(new_gene2.pins))




    return new_gene, new_gene2

def roundifyImage(target):
    radius = image_size/2
    image_center = [image_size/2, image_size/2]
    for i in range(0,target.shape[0]):
        for j in range(0, target.shape[1]):
            dist = sqrt( (i - image_center[0])**2 + (j - image_center[1])**2)
            if(dist >= radius):
                target[i,j] = 255


if (__name__ == "__main__"):
    canvas = np.full((image_size,image_size), 255, np.uint8)
    target = cv2.imread(target_image_name, cv2.IMREAD_GRAYSCALE)

    if(not target_preprocessed):
        target = cv2.resize(target,(image_size,image_size))
        target = cv2.medianBlur(target,7)
        target = cv2.adaptiveThreshold(target,255,cv2.ADAPTIVE_THRESH_MEAN_C,\
                cv2.THRESH_BINARY,9,4)
        roundifyImage(target)
        cv2.imwrite("target_final.png", target)


    pins = generatePins(number_of_pins)
    drawPins(canvas, pins)

    drawLine(canvas, pins[0], pins[1])

    valid_strings = np.zeros((number_of_pins,number_of_pins), np.uint8)
    for i in range(0, number_of_pins):
        for j in range(0, number_of_pins):
            if(i==j):
                valid_strings[i,j] = False
                continue
            temp_canvas = np.full((image_size,image_size), 0, np.uint8)
            pin1 = pins[i]
            pin2 = pins[j]
            point1 = (floor(offset_x + pin1[0]), floor(offset_y + pin1[1]))
            point2 = (floor(offset_x + pin2[0]), floor(offset_y + pin2[1]))
            cv2.line(temp_canvas, point1, point2, 255-line_color, line_thickness)
            my_cost = compareImages(temp_canvas, target)
            #print("i/j " + str(i) +"/" +str(j) + " Cost: " + str(my_cost))
            if(np.count_nonzero(target==0) - my_cost < 2):
                valid_strings[i,j] = False
            else: 
                valid_strings[i,j] = True
                #print("==NOT BANNED  i/j " + str(i) +"/" + str(j))  
    print("Number of valid:     " + str(np.count_nonzero(valid_strings==True)))
    print("Number of non-valid: " + str(np.count_nonzero(valid_strings==False)))
    Gene.valid_strings = valid_strings.copy()

    generation = 1
    #Initialize Genes
    gene_list = [Gene(intial_gene_length) for i in range(0, number_of_genes)]
    print(len(gene_list[-1].pins))
    while(True):
        print("====GeneLength " + str(len(gene_list)) + "==========================================")
        print(generation)
        ###Evaluate genes and normalize values
        print("Evaluating...")
        evaluation = [evaluateGene(gene, target, pins) for gene in gene_list]
        evaluation_abs = evaluation.copy()
        evaluation_sum = sum(evaluation)
        evaluation = [number/evaluation_sum for number in evaluation]
        
        #print(evaluation_abs)
        #print(evaluation)

        ###Sort
        print("Sorting...")
        indeces_sorted = np.argsort(evaluation)
        indeces_sorted = np.flip(indeces_sorted) #Sort in descending order

        #print(indeces_sorted)
        
        ###Printing generation for reference
        best_gene_of_generation = gene_list[indeces_sorted[-1]]
        best_gene_of_generation_score = evaluation[indeces_sorted[-1]]
        worst_gene_of_generation_score = evaluation[indeces_sorted[1]]
        test_length = len(gene_list[indeces_sorted[-1]].pins)

        print("Length of best gene:                 " + str(len(gene_list[indeces_sorted[-1]].pins)))
        print("Best abs score of generation:        " + str(evaluation_abs[indeces_sorted[-1]]))
        print("SecondBest abs score of generation:  " + str(evaluation_abs[indeces_sorted[-2]]))
        print("ThridBest abs score of generation:   " + str(evaluation_abs[indeces_sorted[-3]]))

        for i in range(0, len(gene_list)):
            counter = 0
            js = []
            for j in range(i, len(gene_list)):
                if(gene_list[i].pins == gene_list[j].pins):
                    counter = counter+1
                    js.append(j)
            if(counter>1):
                print("WHAT THE FUCK IT'S THE SAME GENE  " + str(counter-1) + " i/j " + str(i) +"/" +str(js[1:]))


        if(generation % 20 == 0):
            saveGene(best_gene_of_generation, pins, "best", generation, evaluation_abs[indeces_sorted[-1]])
            #saveGene(gene_list[indeces_sorted[-2]], pins, "SecondBest", generation, evaluation[indeces_sorted[-2]])
            #saveGene(gene_list[indeces_sorted[-3]], pins, "ThirdBest", generation, evaluation[indeces_sorted[-3]])

        
        ###Culling
        print("Culling...")
        survival_sum = 0
        i = 0
        while(survival_sum < survival_threshold):
            survival_sum = survival_sum + evaluation[indeces_sorted[i]]
            i = i+1

        
        ###Saving survivors
        survived_genes = [gene_list[indeces_sorted[k]] for k in range(i, len(gene_list))]
        num_of_survivors = len(survived_genes)
        print("Number of survivors  "+ str(num_of_survivors))

        ###Creating new Genes to fill the space
        print("Creating New...")
        while(len(survived_genes) < number_of_genes):
            index1 = random.randint(0, num_of_survivors - 1)
            index2 = random.randint(0, num_of_survivors - 1)
            while(index1 == index2):
                index2 = random.randint(0, num_of_survivors - 1)
            #print("Fusing " + str(index1) + " + " +str(index2))
            new_gene1, new_gene2 = createNewGene(survived_genes[index1], survived_genes[index2])
            survived_genes.append(new_gene1)
            survived_genes.append(new_gene2)

        ###Check to see if mutations/extensions should happen
        print("Applying mutations...")
        mutate = [random.random()<mutate_chance for i in range(0, len(survived_genes))]
        grow = [random.random()<grow_chance for i in range(0, len(survived_genes))]
        [survived_genes[i].mutate() for i in range(0, len(survived_genes)) if mutate[i]]
        [survived_genes[i].grow() for i in range(0, len(survived_genes)) if grow[i]]

        ###Setting up for next round
        gene_list = survived_genes.copy()
        generation = generation + 1
        #break
        


    gene1 = Gene(intial_gene_length)
    gene2 = Gene(intial_gene_length)
    print(gene1.pins)

    print(evaluateGene(gene1, target, pins))
    createNewGene(gene1, gene1)

    drawGene(gene1, pins)

    cv2.imshow("canvas", canvas)
    cv2.imshow("target", target)
    cv2.waitKey(0)

