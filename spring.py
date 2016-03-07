# -*- coding:utf-8 -*-
import re,time,os,math
from os.path import normpath,join
#此文件为3.7日,找雷老师进行讨论后的对原有的C(k)针对碱基进行相关性计算带来的不适用性进行修正的版本
#引入C(d)的方式,计算方法有二:
#首先都是要筛选出距离满足d的所有cpg甲基化值对记为以下格式[[a_1 a_2] [a_1' a_2']......]
#1.采用统计概率的方式统计出a_pre,a_post分别为(0,0):N1,(0,1):N2,(1,0):N3,(1,1):N4四种情况的数量分布,并依照下列公式进行C(d)的计算
#C(d)=(N1/(N1+N2)+N4/(N3+N4))-1,然后画出d与C(d)的相关曲线
#2.采用皮尔逊相关系数r(d)=(sum_{i=0}^n(X_i-\overline{X})*(Y_i-\overline{Y}))/(sqrt((sum_{i=0}^n(X_i-\overline{X})^2)*(sum_{i=0}^n(Y_i-\overline{Y})^2)))
#画出d与r(d)的分布曲线
#从bed文件中读取所有CpG的位置,并记下来,保存为一个hash表,根据key排序,key为int型
#chr_no_list=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X','Y']
chr_no_list=['X','Y']
def read_bed_file_and_store_pos_to_a_struct(bedfile_path):
    struct_to_store={}
    bed_file=open(bedfile_path,'r')
    re_pattern=r'(\d+)\s([\d]+\.[\d]*)\s'
    line=bed_file.readline()
    while line:
        match = re.search(re_pattern,line)
        if match:
            pos=int(match.group(1))
            methy_level=float(match.group(2))
            struct_to_store[pos]=methy_level
        line=bed_file.readline()
    bed_file.close()
    #struct_to_store=sorted(struct_to_store.items(), key=lambda d:d[0])
    return struct_to_store
#根据距离d筛选在hash表中存储的所有满足要求距离的CpG对,输出为[[a_1 a_2] [a_1' a_2']......]形式
def filter_d_length_to_generate_CpG_pairs(CpG_pos_and_methy_struct,d):
    array_to_store_pairs=[]
    for key in CpG_pos_and_methy_struct.keys():
        pos=key
        methy_level_1=CpG_pos_and_methy_struct[key]
        pos_off_by_d=pos+d
        if CpG_pos_and_methy_struct.has_key(pos_off_by_d):
            methy_level_2=CpG_pos_and_methy_struct[pos_off_by_d]
            array_to_store_pairs.append([methy_level_1,methy_level_2])
    return array_to_store_pairs
#根据第一种方法计算C(d)的值
def calc_C_d_by_prob_N1_N2(CpG_pairs,threshold_to_distinct_=0.5):
    N1=N2=N3=N4=0
    for pair in CpG_pairs:
        first=0
        second=0
        if pair[0]>=0.5:
            first=1
        if pair[1]>=0.5:
            second=1
        #N1情况
        if first==0 and second==0:
            N1=N1+1.0
        elif first==0 and second==1:
            N2=N2+1.0
        elif first==1 and second==0:
            N3=N3+1.0
        elif first==1 and second==1:
            N4=N4+1.0
    if (N1+N2)==0 or (N3+N4)==0:
        return -1
    C_d=(N1/(N1+N2)+N4/(N3+N4))-1.0

    return C_d
def calc_C_d_by_pearson_correlation(CpG_pairs):
    sum_pre=0.0
    sum_post=0.0
    for pair in CpG_pairs:
        sum_pre=sum_pre+pair[0]
        sum_post=sum_post+pair[1]
    length=len(CpG_pairs)
    mean_1=sum_pre/length
    mean_2=sum_post/length

    sum_up=0.0
    sum_down_left=0.0
    sum_down_right=0.0
    for pair in CpG_pairs:
        X_i=pair[0]
        Y_i=pair[1]
        sum_up=sum_up+(X_i-mean_1)*(Y_i-mean_2)
        sum_down_left=sum_down_left+(X_i-mean_1)*(X_i-mean_1)
        sum_down_right=sum_down_right+(Y_i-mean_2)*(Y_i-mean_2)
    sum_down=math.sqrt(sum_down_left*sum_down_right)
    if sum_down==0:
        return -1
    r_d=sum_up/sum_down
    return r_d
if __name__ == '__main__':
    sample_name="GSM1386021_2cell_mc_CG_paternal_plus"

    chr_out_dir="chr_out_cd_and_rd"
    if not os.path.exists(chr_out_dir):
        os.makedirs(chr_out_dir)
    for chr_no in chr_no_list:
        bed_file_path=sample_name+os.sep+"chr"+str(chr_no)+".bed"
        CpG_pos_and_methy_struct=read_bed_file_and_store_pos_to_a_struct(bed_file_path)
        d_list=range(2,2000)
        out_C_d_file_name=chr_out_dir+os.sep+"chr"+str(chr_no)+"C_d.csv"
        out_C_d_file=open(out_C_d_file_name,'w')

        out_R_d_file_name=chr_out_dir+os.sep+"chr"+str(chr_no)+"r_d.csv"
        out_R_d_file=open(out_R_d_file_name,'w')
        for d in d_list:
            CpG_pairs=filter_d_length_to_generate_CpG_pairs(CpG_pos_and_methy_struct,d)
            if len(CpG_pairs)==0:
                print "passed d=%d" %d
                continue
            C_d=calc_C_d_by_prob_N1_N2(CpG_pairs)
            if C_d!=-1:
                line=str(d)+","+str(C_d)+"\n"
                out_C_d_file.write(line)
                print "finish chr%s d=%d run , C_d=%f" %(chr_no,d,C_d)

                r_d=calc_C_d_by_pearson_correlation(CpG_pairs)
                if r_d!=-1:
                    line2=str(d)+","+str(r_d)+"\n"
                    out_R_d_file.write(line2)
                    print "finish chr%s d=%d run , r_d=%f" %(chr_no,d,r_d)
        out_C_d_file.close()
        out_R_d_file.close()