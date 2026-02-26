"""
script to parsed a fasta file and create a distance matrix 
""" 
import argparse as ag
import pandas as pd
import numpy as np

def norm_mtrx(mtrx, disimilarity=True): 
    """
    Function that normalize the matrix 
    """
    #get the diogonal of 0 to do not count the slef comparison for the max value 
    for col in mtrx.columns: 
        #we asign the diagonal of the matrix in 0 
        mtrx.loc[col,col]=0
    #first we have to obtain the max value in the matrix 
    max_score = np.nanmax(mtrx.values) 
    #now we normalized the matrix for all the values  
    mtrx_norm= mtrx/max_score 
    
    if disimilarity: 
        mtrx_dis = 1-mtrx_norm   
    else: 
        mtrx_dis= None
    
    mtrx_norm = mtrx_norm.fillna(0)
    mtrx_dis = mtrx_dis.fillna(0)
    return [mtrx_norm,mtrx_dis]


def generate_mtrx(DF_aligns, way):
    """
    Function that generates matrix 
    """ 
    axis_matrix=list(DF_aligns.drop_duplicates(subset="query")["query"]) 

    #we make the matrix using pivot 
    Matrix_df=DF_aligns.pivot(
        index="query",
        columns="subject",
        values="score"
    ) 
    #make sure the axis of the matrix and the comlumns are in the same order 
    Matrix_df=Matrix_df.reindex(index=axis_matrix, columns=axis_matrix) 
    
    #we search what optio does the user want
    match way:
        case "upp": 
            #save a copy of the columns
            columns_mtrx=list(Matrix_df.columns)
            for i in Matrix_df.index:
                for j in columns_mtrx:
                    if pd.notna(Matrix_df.loc[i,j]):
                        Matrix_df.loc[j,i] = Matrix_df.loc[i,j]
                    elif pd.notna(Matrix_df.loc[j,i]):
                        Matrix_df.loc[i,j] = Matrix_df.loc[j,i] 
                #for every column we have changue in order to the index we drop it for do not revert the changues 
                columns_mtrx.remove(i) 
        case "max": 
            #we select the max value compared with all the positions 
            Matrix_df = np.maximum(Matrix_df, Matrix_df.T)
    return Matrix_df

        
    
def parser(): 
    """
    Function to parse the arguments provided by the user
    """ 
    parser=ag.ArgumentParser(description="Parser to store all the args provided by comand line") 
    
    parser.add_argument("-t", "--tsv",
                            default="../results/citrocromeC_seq_hmmr_clean_aling_clean.tsv", 
                            type=str, 
                            required=True, 
                            help="Path of the tsv file with the results of the alignment in the format query /t subject /t score") 
    parser.add_argument("-w", "--way", 
                            type=str, 
                            required=True, 
                            help="Specify the way to build the matrix 'upp' for consider the query result over the subject (upper triangular), 'max' for consider the max result") 
    args = parser.parse_args() 
    return args



def main(): 
    """
    Function that control the general program execution 
    """ 
    
    #first we obtian the args pased by the user 
    arguments=parser()
    #now we can assigned them to a variable 
    flag_way = arguments.way 
    Df_raw = pd.read_csv(arguments.tsv, sep="\t",names=["query", "subject","score"])
    
    #we obtain the matrix 
    mtrx_score=generate_mtrx(Df_raw, flag_way) 
    #now whe normalized 
    mtrx_list=norm_mtrx(mtrx_score, True)
    
    #finally save the files 
    mtrx_normalizada = mtrx_list[0]
    mtrx_distances = mtrx_list[1] 
    mtrx_distances.to_csv("../results/dist_matrix.tsv", sep="\t")
    mtrx_normalizada.to_csv("../results/nrom_matrix.tsv", sep="\t")
    mtrx_score.to_csv("../results/score_matrix.tsv", sep="\t")

if __name__ == "__main__":
	main() 