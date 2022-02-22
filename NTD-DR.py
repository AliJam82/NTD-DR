
import random
from decimal import *
import numpy
numpy.seterr(divide='ignore', invalid='ignore')
import copy
from datetime import datetime
#from openpyxl import load_workbook
from scipy import linalg
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
from sklearn.metrics import accuracy_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import precision_score
#This function reads data and return the input matrices in addition to their sizes
def read_data():
    index = 0
    CT_col_count=0
    CD_row_count=0
    TD_row_count=0
    CT_col_count=0
    CD_col_count=0
    TD_col_count=0
    
    #This block gets the interaction matrix of Drugs and Targets
    with open('Chem_Tar_Association.csv') as f:
        all_lines1 = f.readlines()  # get all the lines
        all_lines1 = [x.strip() for x in all_lines1]
        cols = all_lines1[0].split(',')
        CT_matrix = numpy.zeros((len(all_lines1) - 1, len(cols)-1))
        CT_matrix_T = numpy.transpose(CT_matrix)
        all_lines1 = numpy.delete(all_lines1, 0)
        x_count = 0
        CT_row_count = len(all_lines1)

        for i in range(len(all_lines1)):
            rows1 = all_lines1[i].split(',')
            rows1 = numpy.delete(rows1, 0)
            CT_col_count = len(rows1)
            for j in range(len(rows1)):
                CT_matrix[i][j] = rows1[j]
    x_count = 0


    count_CT = 0
    for i in range(len(all_lines1)):
        for j in range((len(rows1))):
            if CT_matrix[i][j] == 1:
                count_CT += 1




    # This block gets the interaction matrix of Drugs and Diseases
    with open('Chem_Dis_Association.csv') as f:
        all_lines2 = f.readlines()  # get all the lines
        all_lines2 = [x.strip() for x in all_lines2]
        cols = all_lines2[0].split(',')
        CD_matrix = numpy.zeros((len(all_lines2) - 1, len(cols)-1))
        CD_matrix_T = numpy.transpose(CD_matrix)
        all_lines2 = numpy.delete(all_lines2, 0)
        x_count = 0
        CD_row_count = len(all_lines2)

        for i in range(len(all_lines2)):
            rows2 = all_lines2[i].split(',')
            rows2 = numpy.delete(rows2, 0)
            CD_col_count = len(rows2)
            for j in range(len(rows2)):
                CD_matrix[i][j] = rows2[j]
    x_count = 0

    count_CD = 0
    for i in range(len(all_lines2)):
        for j in range(len(rows2)):
            if CD_matrix[i][j] == 1:
                count_CD += 1



    # This block gets the interaction matrix of Targets and Diseases
    with open('Tar_Dis_Association.csv') as f:
        all_lines3 = f.readlines()  # get all the lines
        all_lines3 = [x.strip() for x in all_lines3]
        cols = all_lines3[0].split(',')
        TD_matrix = numpy.zeros((len(all_lines3) - 1, len(cols)-1))
        TD_matrix_T = numpy.transpose(TD_matrix)
        all_lines3 = numpy.delete(all_lines3, 0)
        x_count = 0
        TD_row_count = len(all_lines3)
        for i in range(len(all_lines3)):
            rows3 = all_lines3[i].split(',')
            rows3 = numpy.delete(rows3, 0)
            TD_col_count = len(rows3)
            for j in range(len(rows3)):
                TD_matrix[i][j] = rows3[j]


    count_TD = 0
    for i in range(len(all_lines3)):
        for j in range(len(rows3)):
            if TD_matrix[i][j] == 1:
                count_TD += 1


    chem_len = len(CT_matrix)
    dis_len = len(rows3)
    tar_len = len(TD_matrix)




    # This block gets the similarity matrix of drugs
    with open('Chem_Sim.csv') as f:
        all_lines = f.readlines()
        all_lines = [x.strip() for x in all_lines]

        SC_matrix = numpy.zeros((len(all_lines) - 1, len(all_lines) - 1))
        DC_matrix = numpy.zeros((len(all_lines) - 1, len(all_lines) - 1))
        all_lines = numpy.delete(all_lines, 0)

        for i in range(len(all_lines)):
            rows = all_lines[i].split(',')
            rows = numpy.delete(rows, 0)
            for x in range(len(rows)):
                SC_matrix[i][x] = rows[x]
            index += 1
    x_count = 0

    for x in range(len(all_lines) - 1):
        SC_x = SC_matrix[x]
        for y in range(len(all_lines) - 1):
            x_count += SC_x[y]
        DC_matrix[x][x] = Decimal(x_count)
        x_count = 0

    # This block gets the similarity matrix of targets
    with open('Tar_Sim.csv') as f:
        all_lines = f.readlines()  # get all the lines
        all_lines = [x.strip() for x in all_lines]

        ST_matrix = numpy.zeros((len(all_lines) - 1, len(all_lines) - 1))
        DT_matrix = numpy.zeros((len(all_lines) - 1, len(all_lines) - 1))
        all_lines = numpy.delete(all_lines, 0)

        for i in range(len(all_lines)):
            rows = all_lines[i].split(',')
            rows = numpy.delete(rows, 0)
            for x in range(len(rows)):
                ST_matrix[i][x] = rows[x]
            index += 1
    x_count = 0

    for x in range(len(all_lines) - 1):
        ST_x = ST_matrix[x]
        for y in range(len(all_lines) - 1):
            x_count += ST_x[y]
        DT_matrix[x][x] = Decimal(x_count)
        x_count = 0


    # This block gets the similarity matrix of diseases
    with open('Dis_Sim.csv') as f:
        all_lines = f.readlines()  # get all the lines
        all_lines = [x.strip() for x in all_lines]

        SD_matrix = numpy.zeros((len(all_lines) - 1, len(all_lines) - 1))
        DD_matrix = numpy.zeros((len(all_lines) - 1, len(all_lines) - 1))
        all_lines = numpy.delete(all_lines, 0)

        for i in range(len(all_lines)):
            rows = all_lines[i].split(',')
            rows = numpy.delete(rows, 0)
            for x in range(len(rows)):
                SD_matrix[i][x] = rows[x]
            index += 1
    x_count = 0

    for x in range(len(all_lines) - 1):
        SD_x = SD_matrix[x]
        for y in range(len(all_lines) - 1):
            x_count += SD_x[y]
        DD_matrix[x][x] = Decimal(x_count)
        x_count = 0
    return chem_len,dis_len,tar_len,DC_matrix,SC_matrix,DT_matrix,ST_matrix,DD_matrix,SD_matrix,CT_matrix,CD_matrix,TD_matrix,CT_col_count,CD_col_count,TD_col_count,CT_row_count,CD_row_count,TD_row_count
#------------------------------------------------------------------------
#Your main function that calculates output matrix
def calculation( CT_matrix,CD_matrix,TD_matrix,chem_len,dis_len,tar_len,DC_matrix,SC_matrix,DT_matrix,ST_matrix,DD_matrix,SD_matrix,R, alpha, lambda_CT, lambda_CD, lambda_TD, gamma_C, gamma_T, gamma_D,CT_col_count,CD_col_count,TD_col_count,CT_row_count,CD_row_count,TD_row_count):


    LC = DC_matrix - SC_matrix  #Construct the Laplacian matrix for drugs
    LT = DT_matrix - ST_matrix  #Construct the Laplacian matrix for targets
    LD = DD_matrix - SD_matrix  #Construct the Laplacian matrix for disease


    #This block constructs the drug-target-disease tensor (I drugs * J targets * K disease)
    tensor_X = numpy.zeros((dis_len, chem_len, tar_len)) #This tensor stores triple association information
    tensor_Y = numpy.zeros((dis_len, chem_len, tar_len)) #This tensor stores predicted triplet association scores
    P_matrix_CT = numpy.zeros((chem_len, tar_len)) #This matrix stores predicted drug-target association scores
    P_matrix_CD = numpy.zeros((chem_len, dis_len)) #This matrix stores predicted drug-disease association scores
    P_matrix_TD = numpy.zeros((tar_len, dis_len)) #This matrix stores predicted target-disease association scores


    count_tensor = 0
    for chem in range(chem_len):
        for dis in range(dis_len):
            for tar in range(tar_len):
                if CT_matrix[chem][tar] == CD_matrix[chem][dis] == TD_matrix[tar][dis] == 1:
                    tensor_X[dis][chem][tar] = 1
                    count_tensor += 1


    matrix_array = []
    matrix_array_t = []
    matrix_array_t_flat = []

    for tensor_dis in range(len(tensor_X)):
        matrix_array.append(tensor_X[tensor_dis])
        matrix_array_t.append(numpy.transpose(tensor_X[tensor_dis]))

    X_1 = numpy.hstack(matrix_array)    #Construct the mode-1 matricization of tensor_X
    X_2 = numpy.hstack(matrix_array_t)  # Construct the mode-2 matricization of tensor_X


    for mat in range(len(matrix_array_t)):
        matrix_array_t_flat.append(matrix_array_t[mat].flatten())

    X_3 = numpy.vstack(matrix_array_t_flat)  # Construct the mode-3 matricization of tensor_X


    A = numpy.random.random((CT_row_count, R))
    B = numpy.random.random((CT_col_count, R))
    C = numpy.random.random((TD_col_count, R))
    a_difference_old = numpy.zeros((CT_row_count, R))
    a_difference = numpy.zeros((CT_row_count, R))
    b_difference_old = numpy.zeros((CT_col_count, R))
    b_difference = numpy.zeros((CT_col_count, R))
    c_difference_old = numpy.zeros((TD_col_count, R))
    c_difference = numpy.zeros((TD_col_count, R))

    countera = 0
    counterb = 0
    counterc = 0
    a_cell_count = 0
    b_cell_count = 0
    c_cell_count = 0
    A_flag = True
    B_flag = True
    C_flag = True
    iteratea = 0
    iterateb = 0
    iteratec = 0
    counta = 0
    countb = 0
    countc = 0

    while A_flag:

        # Updating Matrix A
        A_old = numpy.copy(A)

        # a_numerator
        # Term 1
        alpha_X1 = alpha * X_1
        CB = linalg.khatri_rao(C, B)
        alpha_X1_CB = numpy.matmul(alpha_X1, CB)

        # Term 2
        GammaC_SC = gamma_C * SC_matrix
        GammaC_SC_A = numpy.matmul(GammaC_SC, A)

        # Term 3
        lambda_CT_Act = lambda_CT * CT_matrix
        lambda_CT_Act_B = numpy.matmul(lambda_CT_Act, B)

        # Term 4
        lambda_CD_Acd = lambda_CD * CD_matrix
        lambda_CD_Acd_C = numpy.matmul(lambda_CD_Acd, C)

        a_numerator = alpha_X1_CB + GammaC_SC_A + lambda_CT_Act_B + lambda_CD_Acd_C

        # a_denominator
        # Term1
        CT_C = numpy.matmul(numpy.transpose(C), C)
        alpha_CT_C = alpha * CT_C
        BT_B = numpy.matmul(numpy.transpose(B), B)
        alpha_CT_C_BT_B = numpy.multiply(alpha_CT_C, BT_B)

        # Term2
        lambda_CT_BT_B = lambda_CT * BT_B

        # Term3
        lambda_CD_CT_C = lambda_CD * CT_C

        # Term4
        parenthesis_A = alpha_CT_C_BT_B + lambda_CT_BT_B + lambda_CD_CT_C

        # Term 5
        A_parenthesis = numpy.matmul(A, parenthesis_A)

        # Term 5
        gammaC_DC = gamma_C * DC_matrix
        gammaC_DC_A = numpy.matmul(gammaC_DC, A)

        a_denominator = A_parenthesis + gammaC_DC_A

        A = numpy.multiply(A_old, (numpy.divide(a_numerator, a_denominator)))
        a_difference_old = a_difference
        a_difference = numpy.absolute(A - A_old)
        counta += 1




        # The updating matrix A
        for i in range(CT_row_count):
            for j in range(R):
                countera += 1
                if a_difference[i][j] <= 0.0001:
                    a_cell_count += 1
                    if a_cell_count == (CT_row_count * R):
                        iteratea = 1
                        print("Count_a = ", countera)
                    elif i == CT_row_count:
                        a_cell_count = 0
        if iteratea == 1:
            A_flag = False


    while B_flag == True:

        # Updating Matrix B
        B_old = numpy.copy(B)

        # b_numerator
        # Term 1
        alpha_X2 = alpha * X_2
        CA = linalg.khatri_rao(C, A)
        alpha_X2_CA = numpy.matmul(alpha_X2, CA)

        # Term 2
        GammaT_ST = gamma_T * ST_matrix
        GammaT_ST_B = numpy.matmul(GammaT_ST, B)

        # Term 3
        lambda_CT_Act_T = lambda_CT * (numpy.transpose(CT_matrix))
        lambda_CT_Act_T_A = numpy.matmul(lambda_CT_Act_T, A)

        # Term 4
        lambda_TD_Atd = lambda_TD * TD_matrix
        lambda_TD_Atd_C = numpy.matmul(lambda_TD_Atd, C)

        b_numerator = alpha_X2_CA + GammaT_ST_B + lambda_CT_Act_T_A + lambda_TD_Atd_C

        # b_denominator
        # Term 1
        AT_A = numpy.matmul(numpy.transpose(A), A)
        alpha_CT_C_AT_A = numpy.multiply(alpha_CT_C, AT_A)

        # Term 2
        lambda_CT_AT_A = lambda_CT * AT_A

        # Term 3
        lambda_TD_CT_C = lambda_TD * CT_C

        # Term 4
        parenthesis_B = alpha_CT_C_AT_A + lambda_CT_AT_A + lambda_TD_CT_C

        # Term 5
        B_parenthesis = numpy.matmul(B, parenthesis_B)

        # Term 6
        gammaT_DT = gamma_T * DT_matrix
        gammaT_DT_B = numpy.matmul(gammaT_DT, B)

        b_denominator = B_parenthesis + gammaT_DT_B

        B = numpy.multiply(B_old, (numpy.divide(b_numerator, b_denominator)))
        b_difference_old = b_difference
        b_difference = numpy.absolute(B - B_old)
        countb +=1



        # The updating matrix B
        for i in range(CT_col_count):
            for j in range(R):
                counterb += 1
                if b_difference[i][j] <= 0.0001:
                    b_cell_count += 1
                    if b_cell_count == (CT_col_count * R):
                        iterateb = 1
                        print("Count_b = ", counterb)
                    elif i == CT_col_count:
                        b_cell_count = 0
        if iterateb == 1:
            B_flag = False


    while C_flag == True:
        # Updating Matrix C
        C_old = numpy.copy(C)

        # c_numerator
        # Term 1
        alpha_X3 = alpha * X_3
        BA = linalg.khatri_rao(B, A)
        alpha_X3_BA = numpy.matmul(alpha_X3, BA)

        # Term 2
        GammaD_SD = gamma_D * SD_matrix
        GammaD_SD_C = numpy.matmul(GammaD_SD, C)

        # Term 3
        lambda_CD_Acd_T = lambda_CD * (numpy.transpose(CD_matrix))
        lambda_CD_Acd_T_A = numpy.matmul(lambda_CD_Acd_T, A)

        # Term 4
        lambda_TD_Atd_T = lambda_TD * (numpy.transpose(TD_matrix))
        lambda_TD_Atd_T_B = numpy.matmul(lambda_TD_Atd_T, B)

        c_numerator = alpha_X3_BA + GammaD_SD_C + lambda_CD_Acd_T_A + lambda_TD_Atd_T_B

        # c_denominator
        # Term1
        alpha_BT_B = alpha * BT_B
        alpha_BT_B_AT_A = numpy.multiply(alpha_BT_B, AT_A)

        # Term2
        lambda_CD_AT_A = lambda_CD * AT_A

        # Term3
        lambda_TD_BT_B = lambda_TD * BT_B

        # Term4
        parenthesis_C = alpha_BT_B_AT_A + lambda_CD_AT_A + lambda_TD_BT_B

        # Term 5
        C_parenthesis = numpy.matmul(C, parenthesis_C)

        # Term 6
        gammaD_DD = gamma_D * DD_matrix
        gammaD_DD_C = numpy.matmul(gammaD_DD, C)

        c_denominator = C_parenthesis + gammaD_DD_C

        C = numpy.multiply(C_old, (numpy.divide(c_numerator, c_denominator)))
        c_difference_old = c_difference
        c_difference = numpy.absolute(C - C_old)
        countc += 1


        # The updating matrix C
        for i in range(TD_col_count):
            for j in range(R):
                counterc += 1
                if c_difference[i][j] <= 0.0001:
                    c_cell_count += 1
                    if c_cell_count == (TD_col_count * R):
                        iteratec = 1
                        print("Count_c = ", counterc)
                    elif i == TD_col_count:
                        c_cell_count = 0
        if iteratec == 1:
            C_flag = False

    # This block calculates drug-target association scores
    #print ("alpha = ", alpha, ",Lambda_CT = ", lambda_CT, ",Lambda_CD = ", lambda_CD,",Lambda_TD = ", lambda_TD,",Gamma_C = ", gamma_C, ",Gamma_T = ", gamma_T,",Gamma_D = ", gamma_D,"    This is printed on: ", str(datetime.now()))
    for i in range(chem_len):
        for j in range(tar_len):
            xij = 0
            for rr in range(R):
                xij += A[i][rr] * B[j][rr]
            P_matrix_CT[i][j] = xij
    numpy.savetxt('CT_scores(alpha=' + str(alpha) + ',Lambda_CT=' + str(lambda_CT) + ',Lambda_CD=' + str(lambda_CD) + ',Lambda_TD=' + str(lambda_TD) + ',Gamma_C=' + str(gamma_C) + ',Gamma_T=' + str(gamma_T) + ',Gamma_D=' + str(gamma_D) + ').txt', P_matrix_CT, '%0.5f')


    # This block calculates drug-disease association scores
    #print ("alpha = ", alpha, ",Lambda_CT = ", lambda_CT, ",Lambda_CD = ", lambda_CD,",Lambda_TD = ", lambda_TD,",Gamma_C = ", gamma_C, ",Gamma_T = ", gamma_T,",Gamma_D = ", gamma_D,"    This is printed on: ", str(datetime.now()))
    for i in range(chem_len):
        for k in range(dis_len):
            xik = 0
            for rr in range(R):
                xik += A[i][rr] * C[k][rr]
            P_matrix_CD[i][k] = xik
    numpy.savetxt('CD_scores(alpha=' + str(alpha) + ',Lambda_CT=' + str(lambda_CT) + ',Lambda_CD=' + str(lambda_CD) + ',Lambda_TD=' + str(lambda_TD) + ',Gamma_C=' + str(gamma_C) + ',Gamma_T=' + str(gamma_T) + ',Gamma_D=' + str(gamma_D) + ').txt', P_matrix_CD, '%0.5f')


    # This block calculates target-disease association scores
    #print ("alpha = ", alpha, ",Lambda_CT = ", lambda_CT, ",Lambda_CD = ", lambda_CD,",Lambda_TD = ", lambda_TD,",Gamma_C = ", gamma_C, ",Gamma_T = ", gamma_T,",Gamma_D = ", gamma_D,"    This is printed on: ", str(datetime.now()))
    for j in range(tar_len):
        for k in range(dis_len):
            xjk = 0
            for rr in range(R):
                xjk += B[j][rr] * C[k][rr]
            P_matrix_TD[j][k] = xjk
    numpy.savetxt('TD_scores(alpha=' + str(alpha) + ',Lambda_CT=' + str(lambda_CT) + ',Lambda_CD=' + str(lambda_CD) + ',Lambda_TD=' + str(lambda_TD) + ',Gamma_C=' + str(gamma_C) + ',Gamma_T=' + str(gamma_T) + ',Gamma_D=' + str(gamma_D) + ').txt', P_matrix_TD, '%0.5f')


    return P_matrix_CT, P_matrix_CD, P_matrix_TD

#---------------------------------------------------------------------------------
#function for calculate random test indexes based on CV type
def find_test_index(CV, test_version):
    if test_version == 'paiwise':
        fold_num = int(len(index)/CV_num)
        test_index = index[(CV * fold_num):((CV + 1) * fold_num)]
        test_index.sort()
        testPosition = all_position[test_index]
    elif test_version == 'rowwise':
        fold_num = int(len(row_index)/CV_num)
        row_test = row_index[(CV * fold_num):((CV + 1) * fold_num)]
        testPosition = []
        for i in row_test:
            for j in range(0, len(CT_matrix[0])):
                testPosition.append([i, j])
    elif test_version == 'columnwise':
        fold_num = int(len(col_index)/CV_num)
        col_test = col_index[(CV * fold_num):((CV + 1) * fold_num)]
        testPosition = []
        for i in (0, len(CT_matrix)):
            for j in col_test:
                testPosition.append([i, j])
    return numpy.array(testPosition)
def modelEvaluation(real_matrix,predict_matrix,testPosition):
       """
       This function computes the evaluation criteria
       real_matrix
       predict_matrix
       testPosition
       """

       # gathering the test position values in real_matrix and predict_matrix into vectors
       
       real_labels=[]
       predicted_probability=[]
       for i in range(0,len(testPosition)):
           real_labels.append(real_matrix[testPosition[i,0], testPosition[i,1]])
           predicted_probability.append(predict_matrix[testPosition[i,0],testPosition[i,1]])

       real_labels = numpy.array(real_labels)
       predicted_probability=numpy.array(predicted_probability)
       predicted_probability=predicted_probability.reshape(-1,1)
       real_labels=real_labels.reshape(-1,1)
       # computing AUPR criteria
       precision, recall, pr_thresholds = precision_recall_curve(real_labels, predicted_probability)
       aupr_score = auc(recall, precision)

       # computing AUC criteria
       fpr, tpr, auc_thresholds = roc_curve(real_labels, predicted_probability)
       auc_score = auc(fpr, tpr)

       # computing best threshold
       all_F_measure=numpy.zeros(len(pr_thresholds))
       for k in range(0,len(pr_thresholds)):
           if (precision[k]+precision[k])>0:
              all_F_measure[k]=2*precision[k]*recall[k]/(precision[k]+recall[k])
           else:
              all_F_measure[k]=0
       max_index=all_F_measure.argmax()
       threshold =pr_thresholds[max_index]
       #print(threshold)

       # binarize the predited probabilities
       predicted_score=numpy.zeros(len(real_labels))
       predicted_score=numpy.where(predicted_probability > (threshold), 1, 0)

       # computing other criteria
       f=f1_score(real_labels,predicted_score)
       accuracy=accuracy_score(real_labels,predicted_score)
       precision=precision_score(real_labels,predicted_score)
       recall=recall_score(real_labels,predicted_score,predicted_probability)

       # gathering all computed criteria
       results=[auc_score, aupr_score, accuracy, f,precision,recall]
       return results

#---------------------------------------------------------------------
#Reading data
chem_len,dis_len,tar_len,DC_matrix,SC_matrix,DT_matrix,ST_matrix,DD_matrix,SD_matrix,CT_matrix,CD_matrix,TD_matrix,CT_col_count,CD_col_count,TD_col_count,CT_row_count,CD_row_count,TD_row_count=read_data()
#Define variables
seed = 0
CV_num = 10
pos_number = 0
random.seed(seed)
fold_num = 5
all_position = []
index = []
row_index = []
col_index = []
# Need to be specified by user
#Which matrix is your target matrix that you want to apply CV on it. For example CT_matrix. Change it for applying CV on other matrices
CV_in= copy.deepcopy(CT_matrix)
#Percentage for test
split = 0.05
test_version = 'paiwise'
testPosition = 0
if test_version == 'paiwise':
    all_position = []
    for i in range(0, len(CV_in)):
        for j in range(0, len(CV_in[0])):
            pos_number = pos_number + 1
            all_position.append([i, j])
    all_position = numpy.array(all_position)
    index = numpy.arange(0, pos_number)

    zero_indices = random.sample(list(index), int(split*len(index)))
    zeroPosition = all_position[index[zero_indices]]
elif test_version == 'rowwise':
    row_index = numpy.arange(len(CV_in))
    zero_indices = random.sample(list(row_index), int(split*len(row_index)))
    zeroPosition = []
    for i in zero_indices:
        for j in range(0, len(CV_in[0])):
            zeroPosition.append([i, j])
elif test_version == 'columnwise':
    col_index = numpy.arange(len(CV_in[0]))
    zero_indices = random.sample(list(col_index), int(split*len(col_index)))
    zeroPosition = []
    for i in range(0, len(CV_in)):
        for j in zero_indices:
            zeroPosition.append([i, j])

zeroPosition = numpy.array(zeroPosition)
# set the Label values for the test positions to zero
for i in range(0, len(zeroPosition)):
    CV_in[zeroPosition[i, 0], zeroPosition[i, 1]] = 0
#save test positions
numpy.savetxt('BlindPossitions.txt', zeroPosition, delimiter=',')

random.shuffle(index)  # shuffle the indices
random.shuffle(row_index)
random.shuffle(col_index)

#fold_num = (pos_number)// CV_num
r=min([chem_len,dis_len,tar_len])

alphas=[]
lambda_CTs=[]
lambda_CD=1
lambda_TDs=[]
gamma_C=1
gamma_Ts= []
gamma_D=1
bestParam=[]
bestAuC=0
f = open("Tunning.txt", "a")
for alpha in alphas:
        for lambda_CT in lambda_CTs:
            for lambda_TD in lambda_TDs:
                for gamma_T in gamma_Ts:
                    Auc=0
                    Aupr=0
                    acc=0
                    f1=0
                    recall=0
                    precision=0
                    for CV in range(0, CV_num):
                        print('*CV:' + str(CV) + "*\n")
                        # seleting test positions
                        test_index = find_test_index(CV, test_version)

                        test_index.sort()
                        testPosition = all_position[test_index]

                        train_Label= copy.deepcopy(CV_in)

                        # set the Label values for the test positions to zero
                        for i in range(0, len(testPosition)):
                            train_Label[testPosition[i, 0], testPosition[i, 1]] = 0
                        #testPosition = list(testPosition)

                        #I pass train_Label for CT_matrix position, you may pass it for other if you change CV_in
                        P_matrix_CT, P_matrix_CD, P_matrix_TD=calculation(train_Label,CD_matrix,TD_matrix,chem_len,dis_len,tar_len,DC_matrix,SC_matrix,DT_matrix,ST_matrix,DD_matrix,SD_matrix,r,alpha, lambda_CT, lambda_CD, lambda_TD, gamma_C, gamma_T, gamma_D,CT_col_count,CD_col_count,TD_col_count,CT_row_count,CD_row_count,TD_row_count)
                        
                        #print(P_matrix_CT.shape,CV_in.shape )
                        result=modelEvaluation(CV_in,P_matrix_CT,testPosition)
                        Auc=Auc+result[0]
                        Aupr=Aupr+result[1]
                        acc=acc+result[2]
                        f1=f1+result[3]
                        precision=precision+result[4]
                        recall=recall+result[5]
                    Auc=Auc/CV_num
                    Aupr=Aupr/CV_num
                    acc=acc/CV_num
                    f1=f1/CV_num
                    precision=precision/CV_num
                    recall=recall/CV_num
                    print(str(Auc)+','+str(Aupr)+','+str(f1)+','+str(acc)+','+str(precision)+','+str(recall))
                    if Auc>bestAuC:
                        bestAuC=Auc
                        bestParam=[]
                        bestParam.append([alpha,lambda_CT,lambda_TD,gamma_T])
                    f.write('alpha:'+str(alpha)+',lambda_CT:'+str(lambda_CT)+',lambda_TD:'+str(lambda_TD)+',gamma_T:'+str(gamma_T)+',AUC:'+str(Auc)+',AUPR:'+str(Aupr)+',F1:'+str(f1)+',ACC:'+str(acc)+',Precision:'+str(precision)+',Recal:'+str(recall))
                    f.write('\n')
print(bestParam, bestAuC)