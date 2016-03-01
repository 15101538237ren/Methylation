# -*- coding:utf-8 -*-
import re,time,os,math,csv
import shutil
def CG_bp_window_distribution(chr_no,input_file_path,out_file_path,out_interval_geq_path,geq_value,out_interval_leq_path,leq_value,window_size=200):
    #infile:DNA序列文件
    in_file = open(input_file_path,'r')
    out_file = open(out_file_path,'w')
    out_interval_geq=open(out_interval_geq_path,'w')
    out_interval_leq=open(out_interval_leq_path,'w')

    out_header="index,chr,start,end\n"
    out_interval_geq.write(out_header)
    out_interval_leq.write(out_header)

    start_time=time.time()
    line_seq=in_file.readline()
    match = re.search(r'chromosome:([^:]+):(\d*):(\d*):(\d*)',line_seq)
    seq_length=int(match.group(4))
    i=0
    cg_seq=""
    cg_seq_cnt=0
    cg_count_now=0
    one_percent_count=seq_length/100
    first_interval_leq=True
    leq_index=1
    geq_index=1
    bool_cg_seq_enough=False
    while i < seq_length+1:
        if (i % one_percent_count== 0):
            percent = i/one_percent_count
            print "%d percent" % percent


        char_read_in=in_file.read(1)

        if char_read_in=="\n":
            continue

        if cg_seq_cnt < window_size and not bool_cg_seq_enough:
            if cg_seq_cnt!=0:
                now_end=cg_seq[len(cg_seq)-1]
                if char_read_in=="G" and now_end=="C":
                    cg_count_now=cg_count_now+1
                if char_read_in=="C" and now_end=="G":
                    cg_count_now=cg_count_now+1
            cg_seq=cg_seq+char_read_in
            cg_seq_cnt=cg_seq_cnt+1
        #not first 200bp
        else:
            bool_cg_seq_enough=True
            if first_interval_leq and cg_seq[0]!="N":
                cg_start_pos_now=i-window_size+1
                out_interval_leq.write(str(leq_index)+","+chr_no+","+str(cg_start_pos_now)+",")
                leq_index=leq_index+1
                first_interval_leq=False
            cg_start_pos=i-window_size+1
            str_out=str(cg_count_now)+"\n"
            out_file.write(str_out)
            cg_count_pre=cg_count_now
            now_end=cg_seq[len(cg_seq)-1]
            if char_read_in=="G" and now_end=="C":
                cg_count_now=cg_count_now+1
            if char_read_in=="C" and now_end=="G":
                cg_count_now=cg_count_now+1
            if cg_seq[0:2]=="CG" or cg_seq[0:2]=="GC":
                cg_count_now=cg_count_now-1
                #end > geq_value interval
            if cg_count_pre > geq_value and cg_count_now <=geq_value:
                cg_start_pos_now=i-window_size+2
                out_interval_geq.write(str(cg_start_pos_now)+"\n")
                #begin > geq_value interval
            elif cg_count_pre <= geq_value and cg_count_now >geq_value:
                cg_start_pos_now=i-window_size+2
                out_interval_geq.write(str(geq_index)+","+chr_no+","+str(cg_start_pos_now)+",")
                geq_index=geq_index+1
            if cg_count_pre >=leq_value and cg_count_now < leq_value:
                cg_start_pos_now=i-window_size+2
                out_interval_leq.write(str(leq_index)+","+chr_no+","+str(cg_start_pos_now)+",")
                leq_index=leq_index+1
                #begin > geq_value interval
            elif cg_count_pre < leq_value and cg_count_now >=leq_value:
                cg_start_pos_now=i-window_size+2
                out_interval_leq.write(str(cg_start_pos_now)+"\n")
            cg_pre=cg_seq[1:len(cg_seq)]
            cg_seq=cg_pre+char_read_in
        i=i+1
    in_file.close()
    out_file.close()
    out_interval_geq.close()
    out_interval_leq.close()
    print time.time() - start_time
    print seq_length
#输入index,chr,start_site,end_site,筛选出符合filter标准的区间,并按原格式输出到out_file中,count_limit表示需要筛选出的区间个数上限
def filter_interval_and_calc_methy(in_file_path,out_file_path,window_length,filter,seq=",",count_limit=None):
    in_file = open(in_file_path,'r')
    out_file = open(out_file_path,'w')

    start_time=time.time()

    filter_type=filter[0]
    filter_low=0
    filter_high=0

    #bewteen
    if filter_type=="b":
        filter_info=str(filter[1:len(filter)]).split("&")
        filter_low=int(filter_info[0])
        filter_high=int(filter_info[1])
    else:
        filter_low=int(filter[1:len(filter)])
    line=in_file.readline()
    line=in_file.readline()
    pattern="(\d+)"+seq+"(\d+)"+seq+"(\d+)"+seq+"(\d+)\n"
    chr_no=re.search(pattern,line).group(2)
    count_now=0

    if count_limit==None:
        do_not_count=True
    else:
        do_not_count=False
    while line:
        match = re.search(pattern,line)
        if match:
            index=int(match.group(1))
            start_site=int(match.group(3))
            end_site=int(match.group(4))
            length_of_interval=end_site+window_length-start_site

            filter_low_cretiria=filter_type=="<" and length_of_interval < filter_low
            filter_high_cretiria=filter_type==">" and length_of_interval > filter_low
            filter_between_cretiria=filter_type=="b" and length_of_interval > filter_low and length_of_interval < filter_high
            count_limit_cretiria=count_now>=count_limit
            if count_limit_cretiria and not do_not_count:
                break
            else:
                if filter_low_cretiria or filter_high_cretiria or filter_between_cretiria:
                        out_file.write(line)
                        count_now=count_now+1
        line=in_file.readline()
    in_file.close()
    out_file.close()

    print time.time() - start_time

def stat_ck_from_dir(in_dir,out_file_path,sep=","):
    out_file = open(out_file_path,'w')
    k_list=range(0,21)
    k_list2=[25,30,40,50,60,70,80,90,100,150,200,250,300]
    k_list.extend(k_list2)
    line="chr_interval"
    for i in k_list:
        line+=","+str(i)
    line+="\n"
    out_file.write(line)
    filelist = os.listdir(in_dir)
    pattern="([^"+sep+"]+)"+sep+"([^"+sep+"]+)"+sep+"([^"+sep+"]+)"+sep+"([^"+sep+"\n]+)\n"
    for name in filelist:
        if name.endswith(".txt"):
            file_to_read=in_dir +os.sep+ name
            in_file=open(file_to_read,'r')
            line=in_file.readline()
            first=True
            while line:
                match = re.search(pattern,line)
                if match:
                    pos=int(match.group(1))
                    ck=float(match.group(4))
                    #不是最后一个元素
                    if first:
                        line_to_write=name[0:len(name)-4]+sep+str(ck)
                        first=False
                    elif pos!=k_list[len(k_list)-1]:
                        line_to_write=sep+str(ck)
                    else:
                        line_to_write=sep+str(ck)+"\n"
                        first=True
                    out_file.write(line_to_write)
                line=in_file.readline()
            in_file.close()
    out_file.close()
#根据csv输入文件计算各列的均值,除去第一列,第一行的选项为first_column_useless,header
def calc_mean_from_csv(input_file_path,header=True,first_column_useless=True,sep=","):
    input_file=open(input_file_path,'r')
    if header==True:
        input_file.readline()
    #保存各列的均值list
    means=[]

    line=input_file.readline()
    first_line=True
    line_count=0
    while line:
        str_list=line.split(sep)
        for index,val in enumerate(str_list):
            #第一行要append
            if index==0:
                if first_column_useless:
                    continue

            float_val=float(val)

            if first_line:
                means.append(float_val)
            else:
                if first_column_useless:
                    means[index-1]=means[index-1]+float_val
                else:
                    means[index]=means[index]+float_val
        if first_line:
            first_line=False
        line=input_file.readline()
        line_count=line_count+1.0
    input_file.close()
    for index,val in enumerate(means):
        means[index]=val/line_count
    input_file_name_str_list=input_file_path.split(".")
    src=input_file_path
    dst=input_file_name_str_list[0]+"_mean."+input_file_name_str_list[1]
    if os.path.exists(dst):
        os.remove(dst)
    shutil.copyfile(src,dst)
    out_file=open(dst,'a')
    out_file.write("mean")
    for index,val in enumerate(means):
        out_file.write(","+str(means[index]))
    out_file.write("\n")
    out_file.close()
    return means
if __name__ == '__main__':
    data_dir="data"
    seq_file_name="Mus_musculus.GRCm38.dna.chromosome.1.fa"

    input_file_path=data_dir+os.sep+seq_file_name
    out_dir="out"
    out_file_tags=["cg_count_40_plus","cg_count_5_minus"]
    for out_file_tag in out_file_tags:
        out_output=out_dir+os.sep+out_file_tag
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        file=out_dir+os.sep+"ck_stat_"+out_file_tag+".csv"
        stat_ck_from_dir(out_output,file,sep=",")
        means=calc_mean_from_csv(file,True,True,",")
    print "Done handle means append!"
        #file_input=open(file,'r')
        #统计输入文件的均值,2-最后一行,2-最后一列的均值.




    # out_file_name="test.csv"
    # out_file_path=out_dir+os.sep+out_file_name
    # out_interval_geq_name="interval_40_plus.csv"
    # out_interval_leq_name="interval_5_minus.csv"
    # out_interval_geq_path=out_dir+os.sep+out_interval_geq_name
    # out_interval_leq_path=out_dir+os.sep+out_interval_leq_name
    # interval_file_name="out_interval_final.csv"
    # out_interval_final_path=out_dir+os.sep+interval_file_name
    # chr_no="1"
    # #CG_bp_window_distribution(chr_no,input_file_path,out_file_path,out_interval_geq_path,40,out_interval_leq_path,5,200)
    #
    # filter_interval_and_calc_methy(out_interval_leq_path,out_interval_final_path,200,">200",",")




