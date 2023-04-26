import pyVCT
import wandb
import skimage
from skimage.measure import label, regionprops, regionprops_table

import feret
from feret.main import Calculater

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from tqdm.notebook import tqdm
from utils.draw import make_image
from numpngw import write_png
from scipy.stats import entropy


def convert_area(exp_area):
    """
    This function just brings the experimental data to the same dimension, which gives the model in the output
    """
    return exp_area * ((50 / 180 * 0.001)**2 / 0.0025**2)


def make_labeled(c_tags, types):
    """
    :params c_tags, types: first two result of pyVCT.py_cpmfem
    :returns: twp np.array with pixels labeled with the cell types 
    """
    CM_labeled=np.zeros((len(c_tags),len(c_tags[0])), dtype=int)
    FB_labeled=np.zeros((len(c_tags),len(c_tags[0])), dtype=int)
    for i in range(len(c_tags)):
        for j in range(len(c_tags[0])):
            if c_tags[i][j]!=0:
                if types[c_tags[i][j]-1]==1:
                    CM_labeled[i][j]=int(c_tags[i][j])
                if types[c_tags[i][j]-1]==2:
                    FB_labeled[i][j]=int(c_tags[i][j])
    return CM_labeled, FB_labeled


def make_df(c_tags, types):
    """
    :params c_tags, types: first two result of pyVCT.py_cpmfem
    
    :returns: returns two pd.Dataframe with data about cm / fb from the model
    """
    CM_labeled, FB_labeled = make_labeled(c_tags, types)
    labels = set(np.unique(CM_labeled)) - {0}
    length=np.zeros(len(labels))
    width=np.zeros(len(labels))
    hull_number=np.zeros(len(labels))
    for i in enumerate(labels):
        feret_calc = Calculater(CM_labeled==i[1], edge=True)
        feret_calc.find_convexhull()
        hull_number[i[0]]=len((feret_calc.hull).T)
        length[i[0]]=feret.max(CM_labeled==i[1], edge=True)
        width[i[0]]=feret.max90(CM_labeled==i[1], edge=True)
    elongation=length/width
    cell_table_CM= regionprops_table(CM_labeled, properties=['label', 'area', 'convex_area'])
    cell_df_CM = pd.DataFrame(cell_table_CM)
    cell_df_CM['length']=length
    cell_df_CM['width']=width
    cell_df_CM['elongation']=elongation
    cell_df_CM['hull_number']=hull_number
    
    labels = set(np.unique(FB_labeled)) - {0}
    length=np.zeros(len(labels))
    width=np.zeros(len(labels))
    hull_number=np.zeros(len(labels))
    for i in enumerate(labels):
        feret_calc = Calculater(FB_labeled==i[1], edge=True)
        feret_calc.find_convexhull()
        hull_number[i[0]]=len((feret_calc.hull).T)
        length[i[0]]=feret.max(FB_labeled==i[1], edge=True)
        width[i[0]]=feret.max90(FB_labeled==i[1], edge=True)
    elongation=length/width
    cell_table_FB= regionprops_table(FB_labeled, 
                            properties=['label', 'area', 'convex_area'])
    cell_df_FB = pd.DataFrame(cell_table_FB)
    cell_df_FB['length']=length
    cell_df_FB['width']=width
    cell_df_FB['elongation']=elongation
    cell_df_FB['hull_number']=hull_number
    
    return cell_df_CM, cell_df_FB


def boulder(img):
    """
    :param img: model simulation picture, result of make_image function
    
    :returns: part_getero - percent of geterogenous cells in the image
    """
    GR_FB=set(np.array([127, 127, 127]))
    GR_CM=set(np.array([127,   0,   0]))
    FB=set(np.array([255, 255, 255]))
    CM=set(np.array([255,   0,   0]))
    MD=set(np.array([0,   0,   0]))
    FIB=set(np.array([ 25,  25,  25]))
    len_gomo_CM=0
    len_gomo_FB=0
    len_getero=0
#count=[GR_FB,GR_CM,FB,CM,MD]
    count=[0,0,0,0,0]
    n=img.shape[0]
    for i in range(1,n-1):
        for j in range(1,n-1):
            neigh=[[i-1,j-1],[i,j-1],[i+1,j-1],[i-1,j],[i+1,j],[i-1,j+1],[i,j+1],[i+1,j+1]]
        
            for k in range(8):
                neighbour = set(img[neigh[k][0]][neigh[k][1]])
                if neighbour==GR_FB:
                    count[0]+=1
                elif neighbour==GR_CM:
                    count[1]+=1
                elif neighbour==FB:
                    count[2]+=1
                elif neighbour==CM:
                    count[3]+=1
                elif neighbour==MD or neighbour==FIB:
                    count[4]+=1
            if set(img[i][j])==GR_CM:
                if count[0]!=0:
                    len_getero+=0.5
                elif count[4]<=1:
                    if count[3] != 0:
                        if count[1]>2:
                            len_gomo_CM+=0.5
                        if count[1]==2:
                            len_gomo_CM+=1
                    
           #dopisat 
            if set(img[i][j])==GR_FB:
                if count[1]!=0:
                    len_getero+=0.5
                elif count[4]<=1:
                    if count[2] != 0:
                        if count[0]>2:
                            len_gomo_FB+=0.5
                        if count[0]==2:
                            len_gomo_FB+=1
            #dopisat
            count=[0,0,0,0,0]
             
        


    part_CM=len_gomo_CM/(len_gomo_CM+len_gomo_FB+len_getero)
    part_FB=len_gomo_FB/(len_gomo_CM+len_gomo_FB+len_getero)
    part_getero=len_getero/(len_gomo_CM+len_gomo_FB+len_getero)
    return part_CM, part_FB, part_getero
  
        
def get_exp_distributions(all_data, cell_day_type):
    """
    :param all_data: pd.DataFrame with the next columns: 
        cell_type - type of the cell(lc = layer cm, lfb = layer fb),
        smooth_area - area of the cell in the model dimensionality
        area_quotient - quotient of smooth area to the rectangle area of the cell
        elongation - quotient of length and width of the cell, always less then 1
        n_podium - number of the cell podiums
    :param cell_day_type: str, p2-1 or d20
    
    :returns: two dicts - distribution of provided parameters in cardiomyocytes and fibroblasts
    """
    if cell_day_type not in all_data.cell_day_type.unique():
        print('Such cells not provided in the data')
        
    data = all_data[all_data.cell_day_type == cell_day_type]

    cm_data = {
        'hull_number_distribution': (data
                                        .where(data.cell_type == 'lc')
                                        .n_podium
                                        .dropna()
                                        .to_numpy()
                                     ),
        'area_quotient_distribution': (data
                                        .where(data.cell_type == 'lc')
                                        .area_quotient
                                        .dropna()
                                        .to_numpy()
                                       ),
        'elongation_distribution': (data
                                        .where(data.cell_type == 'lc')
                                        .elongation
                                        .dropna()
                                        .to_numpy()
                                       ),
        'area_distribution': (data
                                        .where(data.cell_type == 'lc')
                                        .smooth_area
                                        .dropna()
                                        .to_numpy()
                                       )
    }
    

    fb_data = {
        'hull_number_distribution': (data
                                     .where(data.cell_type == 'lfb')
                                     .n_podium
                                     .dropna()
                                     .to_numpy()
                                     ),
        'area_quotient_distribution': (data
                                        .where(data.cell_type == 'lfb')
                                        .area_quotient
                                        .dropna()
                                        .to_numpy()
                                       ),
        'elongation_distribution': (data
                                        .where(data.cell_type == 'lfb')
                                        .elongation
                                        .dropna()
                                        .to_numpy()
                                       ),
        'area_distribution': (data
                                        .where(data.cell_type == 'lfb')
                                        .smooth_area
                                        .dropna()
                                        .to_numpy()
                                       )
    }
    
    return cm_data, fb_data
       

def compute_metrics(cm_model, fb_model, ctags, types, fibers, contacts, cm_exp, fb_exp,
                    cm_elongation_r=(0, 1), cm_area_quotient_r=(0, 1), cm_n_podium_r=(5, 25), cm_area_r=(50, 400),
                    fb_elongation_r=(0, 1), fb_area_quotient_r=(0, 1), fb_n_podium_r=(6, 17), fb_area_r=(0, 400),
                   ):
    """
    :params cm / fb_model: output of make_df function
    :ctags, types, fibers, contacts: output of pyVCT.py_cpmfem
    :cm_exp, fb_exp: experemental distributions, output of get_exp_distributions funtion
    
    :returns: dict with the observed metrics
    """
    eps = 1e-3
    metrics = {}

    img = make_image(ctags, types, fibers, contacts, 0)

    metrics['cm_elongation_entropy'] = entropy(
        np.histogram(cm_exp['elongation_distribution'], density=True, range=cm_elongation_r)[0]+eps,
        np.histogram(1 / cm_model.elongation, density=True, range=cm_elongation_r)[0]+eps
    )
    metrics['cm_area_quotient_entropy'] = entropy(
        np.histogram(cm_exp['area_quotient_distribution'], density=True, range=cm_area_quotient_r)[0]+eps,
        np.histogram(cm_model.area / cm_model.convex_area, density=True, range=cm_area_quotient_r)[0]+eps
    )
    metrics['cm_hull_number_entropy'] = entropy(
        np.histogram(cm_exp['hull_number_distribution'], density=True, range=cm_n_podium_r)[0]+eps,
        np.histogram(cm_model.hull_number, density=True, range=cm_n_podium_r)[0]+eps
    )
    metrics['cm_area_entropy'] = entropy(
        np.histogram(np.apply_along_axis( lambda x: convert_area(x), arr=cm_exp['area_distribution'], axis=0), density=True, range=cm_area_r)[0]+eps,
        np.histogram(cm_model.area, density=True, range=cm_area_r)[0]+eps
    )

    metrics['fb_elongation_entropy'] = entropy(
        np.histogram(fb_exp['elongation_distribution'], density=True, range=fb_elongation_r)[0]+eps,
        np.histogram(1 / fb_model.elongation, density=True, range=fb_elongation_r)[0]+eps
    )
    metrics['fb_area_quotient_entropy'] = entropy(
        np.histogram(fb_exp['area_quotient_distribution'], density=True, range=fb_area_quotient_r)[0]+eps,
        np.histogram(fb_model.area / fb_model.convex_area, density=True, range=fb_area_quotient_r)[0]+eps
    )
    metrics['fb_hull_number_entropy'] = entropy(
        np.histogram(fb_exp['hull_number_distribution'], density=True, range=fb_n_podium_r)[0]+eps,
        np.histogram(fb_model.hull_number, density=True, range=fb_n_podium_r)[0]+eps
    )
    metrics['fb_area_entropy'] = entropy(
        np.histogram(np.apply_along_axis( lambda x: convert_area(x), arr=fb_exp['area_distribution'], axis=0), density=True, range=fb_area_r)[0]+eps,
        np.histogram(fb_model.area, density=True, range=fb_area_r)[0]+eps
    )


    metrics['cm_mean_entropy'] = np.mean(
        [
            metrics['cm_hull_number_entropy'],
            metrics['cm_elongation_entropy'],
            metrics['cm_area_quotient_entropy'],
            metrics['cm_area_entropy'],
        ]
    )
    metrics['fb_mean_entropy'] = np.mean(
        [
            metrics['fb_hull_number_entropy'],
            metrics['fb_elongation_entropy'],
            metrics['fb_area_quotient_entropy'],
            metrics['fb_area_entropy']
        ]
    )

    metrics['mean_entropy'] = np.mean(
        [
            metrics['cm_mean_entropy'],
            metrics['fb_mean_entropy']
        ]
    )
    
    part_CM, part_FB, part_getero = boulder(img)
    if (part_getero < 0.1) or (part_CM < 0.1) or (part_FB < 0.1):
        # for the cases whe simultion is non-physical
        metrics['mean_entropy'] = 100
    metrics['part_getero'] = part_getero
    metrics['part_CM'] = part_CM
    metrics['part_FB'] = part_FB
    
    metrics['image'] = wandb.Image(img)
    
    return metrics


def compare_distributions(cm_model, fb_model, cm_exp, fb_exp):
    """
    :param cm / fb_model: data from model, result of make_df function
    :param cm_exp / fb_exp: data from experiment, result of get_exp_distributions function
    
    plots compration of the exp and model distributions
    """
    names = ['CM', 'FB']
    i = 0
    for model_data, exp_data in [[CM, cm_exp], [fb_model, fb_exp]]:
        plt.figure(figsize=(15,5))
        
        ax = plt.subplot(1, 3, 1)
        sns.histplot(
            1 / (model_data.area / model_data.convex_area),
            label='model',
            kde=True,
            stat='probability'
        );
        sns.histplot(
            1 / exp_data['area_quotient_distribution'],
            label = 'experiment',
            stat='probability',
            kde=True
        );
        plt.xlabel('area')
        plt.legend()
        plt.title(f'{names[i]}')

                                 
        ax = plt.subplot(1, 3, 2)
        sns.histplot(
            model_data.elongation,
            label='model',
            kde=True,
            stat='probability'
        );
        sns.histplot(
            exp_data['elongation_distribution'],
            label = 'experiment',
            kde=True,
            stat='probability'
        ); 
        plt.legend()
        plt.title(f'{names[i]}')
                                 
        ax = plt.subplot(1, 3, 3)
        sns.histplot(
            model_data.hull_number,
            label='model',
            kde=True,
            stat='probability'
        );
        sns.histplot(
            exp_data['hull_number_distribution'],
            label = 'experiment',
            kde=True,
            stat='probability'
        ); 
        plt.legend()
        plt.title(f'{names[i]}')
        plt.show()
        
        i += 1
