# -*- coding:utf-8 -*-
#输入的是序列文件,以及甲基化文件,输出的是序列对甲基化的文件(删除Gap之后)
#跳过所有N的序列同时mapping上去,提取一行关键字,bed 前编号+1=正式maping 的编号,并计算出平均的甲基化值
import re,time,os,math,shutil
from os.path import normpath,join
#chr_no_list=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X','Y']
chr_no_list=['1']
def scan(input_seq_path,input_methy_path,output_file_path,statistics_file,interval_left,interval_right,sep=","):
    seq_file = open(input_seq_path,'r')
    methy_file = open(input_methy_path,'r')
    out_file = open(output_file_path,'w')
    stat_file=open(statistics_file,'w')
    line = methy_file.readline()
    methy_patern=r'(\d+)(\s)([\d]+\.[\d]*)'
    match = re.search(methy_patern,line)
    methy_pos=int(match.group(1))
    value=float(match.group(3))

    pos=0
    char_readin=""
    out_str=""
    done=False

    all_base_count=0
    all_base_value=0.0
    methy_cnt=0

    #remove_the_first_line
    start=time.clock()
    line_seq=seq_file.readline()
    match = re.search(r'chromosome:([^:]+):(\d*):(\d*):(\d*)',line_seq)
    length=int(math.ceil(int(match.group(4))/60.0+1.0))
    i=1
    percent_100=length/100
    while methy_pos<interval_left:
        line = methy_file.readline()
        if line:
            match = re.search(methy_patern,line)
            methy_pos=int(match.group(1))
            value=float(match.group(3))
    while i < length+1:
        line_seq = seq_file.readline()
        # if (i % percent_100== 0):
        #     percent = i/percent_100
        #     print "%d percent" % percent

        line_len=len(line_seq)
        pre_pos=(i-1)*60
        now_pos=(i-1)*60+line_len-1
        if interval_left and interval_right:
            if interval_right-interval_left<60:
                print "interval length must greater than 60!"
                break
            else:
                bool_left_out=now_pos<interval_left
                bool_left_part_in=now_pos>=interval_left and pre_pos < interval_left
                bool_in_or_part_right=pre_pos >= interval_left and pre_pos<interval_right
                booo_right_out=pre_pos>=interval_right
                if bool_left_out:
                    i=i+1
                elif bool_left_part_in or bool_in_or_part_right:
                    for j in range(1,line_len):
                        pos =(i-1)*60+j
                        if pos >=interval_left and pos<=interval_right:
                            char_readin = line_seq[j-1]
                            # if pos==195371430:
                            #     print "haha,195371429"
                            if (pos==methy_pos+1) and (not done):
                                methy_cnt = methy_cnt+1
                                out_str=str(pos)+sep+char_readin+sep+str(value)+"\n"
                                #print "%d %f" %(methy_pos+1,value)
                                out_file.write(out_str)
                                all_base_count=all_base_count+1
                                all_base_value=all_base_value+value
                                line = methy_file.readline()
                                if line:
                                    match = re.search(r'(\d+)(\s)([\d]+\.[\d]*)',line)
                                    methy_pos=int(match.group(1))
                                    value=float(match.group(3))
                                else:
                                    done=True
                            elif char_readin !='N':
                                out_str=str(pos)+sep+char_readin+sep+str(0.0)+"\n"
                                out_file.write(out_str)
                                all_base_count=all_base_count+1
                    i=i+1
                elif booo_right_out:
                    break
        else:
            break
    seq_file.close()
    methy_file.close()
    out_file.close()
    average_value=all_base_value/all_base_count
    finish=time.clock()
    # print "Sum basecount:%d,methy_cnt %d,sum_value %f,average %f,time:%f" % (all_base_count,methy_cnt,all_base_value,average_value,(finish-start))
    stat_file.write(str(all_base_value)+"\n")
    stat_file.write(str(all_base_count)+"\n")
    stat_file.write(str(average_value)+"\n")
    return all_base_value,all_base_count,average_value
#
def calculate(k,a_avg,n,stat_file1,stat_file2,out_file,in_seq=",",sep=","):
    sum_up=0.0
    sum_down=0.0

    reader=stat_file1
    k_reader=stat_file2

    start=time.clock()
    per_10=n/100
    for i1 in range(1,k+1):
        k_reader.readline()

    i=1
    while i< n-k+1:
        # if (i % per_10 == 0):
        #     percent = i/per_10
        #     print "%d step: %d percent" % (k,percent)
        line_k=k_reader.readline()
        line=reader.readline()
        a_ipk=float(line_k.split(in_seq)[2])
        a_i=float(line.split(in_seq)[2])
        sum_up=sum_up+(a_i-a_avg)*(a_ipk-a_avg)
        sum_down=sum_down+(a_i-a_avg)*(a_i-a_avg)
        i=i+1
    while i<n+1:
        # if (i % per_10 == 0):
        #     percent = i/per_10
        #     print "%d step: %d percent" % (k,percent)
        line=reader.readline()
        if line:
            a_i=float(line.split(in_seq)[2])
            sum_down=sum_down+(a_i-a_avg)*(a_i-a_avg)
        i=i+1
    if sum_down!=0.0:
        C_k=sum_up/sum_down
    else:
        C_k=0.0
    str_line=str(k)+sep+str(sum_up)+sep+str(sum_down)+sep+str(C_k)+"\n"
    out_file.write(str_line)
    finish=time.clock()
    print "finish in %f" % (finish-start)
def bed_data_extract_to_methy(chr_no,in_file_path,outfile_path):
    raw_file = open(in_file_path,'r')
    methy_file = open(outfile_path,'w')
    chr_no=str(chr_no)

    line=raw_file.readline()
    count=0
    find=0
    while line:
        count=count+1
        if count % 10000==0:
            print "%d of %s was processed" %(count,in_file_path)
        match = re.search(r'chr'+chr_no+r'\t(\d+)(\s)(\d+)(\s)([\d]+\.[\d]*)(\s)(\d+)(\s)(\d+)',line)
        if match and find==0:
            methy_pos=int(match.group(1))
            value=float(match.group(5))
            out_str=str(methy_pos)+"\t"+str(value)+"\n"
            methy_file.write(out_str)
            find=1
        elif match and find==1:
            methy_pos=int(match.group(1))
            value=float(match.group(5))
            out_str=str(methy_pos)+"\t"+str(value)+"\n"
            methy_file.write(out_str)
        elif find==1 and not match:
            break
        line=raw_file.readline()
    methy_file.close()
    raw_file.close()
    print "finish %s chr%s data processing!" %(in_file_path,chr_no)
def batch_bed_to_methy(in_file_path,out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for chr_no in chr_no_list:
        bed_data_extract_to_methy(chr_no,in_file_path,out_dir+"/chr"+chr_no+".bed")
def extract_position(in_file_path,out_file_path):
    in_file = open(in_file_path,'r')
    in_file2 = open(in_file_path,'r')
    out_file = open(out_file_path,'w')
    pattern=r'(\d+)\t([\d]+\.[\d]*)(\s)'
    line=in_file.readline()
    in_file2.readline()
    match = re.search(pattern,line)
    pos1=int(match.group(1))
    pos2=0
    count=0
    while line:
        line=in_file.readline()
        in_file2.readline()
        count=count+1

        if count%10000==0:
            print "%s have read %d lines" % (in_file_path,count)
        match = re.search(pattern,line)
        pos2=int(match.group(1))

        distance=str(pos2-pos1)
        out_file.write(distance+"\n")

        pos1=pos2
        if not in_file2.readline():
            break
    out_file.close()
    in_file.close()
#查看距离的分布情况是否满足filter输出distance以及两个Cpg位点的位置信息
def calc_distance_filter_and_summary(in_file_path,out_file_path,filter_info,chr_no,pattern,sep=","):
    in_file = open(in_file_path,'r')
    out_file = open(out_file_path,'w')

    line=in_file.readline()
    filter_type=filter_info[0]
    filter_distance=0
    filter_distance2=0

    #bewteen
    if filter_type=="b":
        filter_info=str(filter_info[1:len(filter_info)]).split("&")
        filter_distance=int(filter_info[0])
        filter_distance2=int(filter_info[1])
    else:
        filter_distance=int(filter_info[1:len(filter_info)])

    match = re.search(pattern,line)
    pos1=0
    pos2=int(match.group(1))
    count=0
    id=0

    while line:
        if pos1!=0:
            distance=pos2-pos1

            if filter_type=='<':
                if distance<filter_distance:
                    id=id+1
                    out_file.write(str(id)+sep+str(distance)+sep+chr_no+sep+str(pos1)+sep+str(pos1+1)+sep+str(pos2)+sep+str(pos2+1)+"\n")
            elif filter_type =='>':
                 if distance>filter_distance:
                    id=id+1
                    out_file.write(str(id)+sep+str(distance)+sep+chr_no+sep+str(pos1)+sep+str(pos1+1)+sep+str(pos2)+sep+str(pos2+1)+"\n")
            elif filter_type=='=':
                if distance==filter_distance:
                    id=id+1
                    out_file.write(str(id)+sep+str(distance)+sep+chr_no+sep+str(pos1)+sep+str(pos1+1)+sep+str(pos2)+sep+str(pos2+1)+"\n")
            elif filter_type=='b':
                if distance>=filter_distance and distance<=filter_distance2:
                    id=id+1
                    out_file.write(str(id)+sep+str(distance)+sep+chr_no+sep+str(pos1)+sep+str(pos1+1)+sep+str(pos2)+sep+str(pos2+1)+"\n")
        line=in_file.readline()
        if line!='':
            match = re.search(pattern,line)
            pos1=pos2
            pos2=int(match.group(1))
        if count%10000==0:
            print "%s have calculate %d lines" % (in_file_path,count)
        count=count+1
    out_file.close()
    in_file.close()
def distinct_dis_info(in_file_path,out_file_path,sep=","):
    in_file = open(in_file_path,'r')
    out_file = open(out_file_path,'w')
    line=in_file.readline()
    pattern=r'(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)'

    match=re.search(pattern,line)
    chr_no=match.group(3)
    pos_set=set()
    count=0
    while line:
        count=count+1
        if count%1000==0:
            print "%s have disticted %d lines" % (in_file_path,count)
        match=re.search(pattern,line)
        pos1=int(match.group(4))
        pos2=int(match.group(6))
        pos_set.add(pos1)
        pos_set.add(pos2)
        line=in_file.readline()
    pos_list=list(pos_set)
    pos_list.sort()
    index=0
    out_file.write("index,chr_no,start,end\n")
    for item in pos_list:
        index=index+1
        out_file.write("site"+str(index)+sep+chr_no+sep+str(item)+sep+str(item+1)+"\n")
    out_file.close()
    in_file.close()
def calc_distance_filter_and_summary_from_dir(in_file_dir,out_file_dir,filter_info,chr_no,pattern,append=".dat"):
    if not os.path.exists(out_file_dir):
        os.makedirs(out_file_dir)
    pattern=r'(\d+)\t([\d]+\.[\d]*)(\s)'
    for chr_no in chr_no_list:
        calc_distance_filter_and_summary(in_file_dir+os.sep+"chr"+chr_no+".bed",out_file_dir+os.sep+"chr"+chr_no+"_dis_info"+append,filter_info,chr_no,pattern)
#把文件夹中,chr_no_list的bed文件里面所有的Position之间的距离,分开输出到out_dir文件夹
def extract_position_from_dir(in_dir,out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for chr_no in chr_no_list:
        extract_position(in_dir+os.sep+"chr"+chr_no+".bed",out_dir+os.sep+"chr"+chr_no+"_dis.csv")

def count_methylation_distribute(in_file_path,out_file_path):
    in_file = open(in_file_path,'r')
    out_file = open(out_file_path,'w')

    print "start to count_methylation_distribute of %s" % in_file_path
    dict_count={}
    line=in_file.readline()
    pattern=r'(\d+)(\s)'
    while line:
        match=re.match(pattern,line)
        num_count_now=int(match.group(1))

        if dict_count.has_key(num_count_now):
            dict_count[num_count_now]=int(dict_count.get(num_count_now))+1
        else:
            dict_count[num_count_now]=1
        line=in_file.readline()
    dict_out=sorted(dict_count.items(), key=lambda d:d[0])
    out_file.write("distance,count\n")
    for tuple_item in dict_out:
        # if tuple_item[0]<2000:
        #
        # else:
        #     break
        out_file.write(str(tuple_item[0])+","+str(tuple_item[1])+"\n")
    out_file.close()
    in_file.close()
    print "finished count_methylation_distribute of %s" % in_file_path
#计算文件夹内所有chr_no_list的距离的分布情况,输出到out_dir
def count_methylation_distribute_from_dir(in_dir,out_dir,in_append=".csv",out_append=".dat"):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for chr_no in chr_no_list:
        count_methylation_distribute(in_dir+os.sep+"chr"+chr_no+"_dis"+in_append,out_dir+os.sep+"chr"+chr_no+"_count"+out_append)
def merge_count_distribute_from_dir(in_dir,out_file_path,in_append=".csv"):
    out_file = open(out_file_path,'w')
    dict_count={}
    for chr_no in chr_no_list:
        infile=open(in_dir+os.sep+"chr"+chr_no+"_count"+in_append,'r')
        print "start to merge %s" % infile.name
        line=infile.readline()
        pattern=r'(\d+),(\d+)(\s)'
        while line:
            match=re.match(pattern,line)
            key_match=int(match.group(1))
            value_match=int(match.group(2))
            if dict_count.has_key(key_match):
                dict_count[key_match]=dict_count.get(key_match)+value_match
            else:
                dict_count[key_match]=value_match
            line=infile.readline()
        print "finished to merge %s" % infile.name
        infile.close()
    dict_out=sorted(dict_count.items(), key=lambda d:d[0])
    for items in dict_out:
        # if items[0]<=2000:
        #
        # else:
        #     break
        out_file.write(str(items[0])+","+str(items[1])+"\n")
    out_file.close()
#获得所有基因的启动子的位置(根据输入文件,上下游的base数)
def get_all_gene_promoters(in_file_path,out_file_path,upstream,downstream,sep=","):
    in_file = open(in_file_path,'r')
    out_file = open(out_file_path,'w')

    line=in_file.readline()
    print "start to get promoters of %s" % in_file.name
    pattern=r'(\d+|X|Y)\s([^\s]+)\s(gene)\s(\d+)\s(\d+)\s([^\s]+)\s([\+|\-])\s([^\s]+)\s([^\s]+)\s([^\s]+)\s([^\s]+)\s([^\s]+)\s([^\s]+)\s([^\s]+)'
    count=0
    while line:
        match=re.match(pattern,line)
        if match:
            chr_no=match.group(1)
            TSS=int(match.group(4))
            #TTS=int(match.group(8))
            promoter_S=TSS-upstream
            promoter_T=TSS+downstream
            strand=match.group(7)
            gene_name=match.group(14)
            gene_name=gene_name[1:len(gene_name)-2]
            out_file.write(chr_no+sep+str(promoter_S)+sep+str(promoter_T)+sep+strand+sep+gene_name+"\n")
        if count%10000==0:
            print "%s have read %d lines" % (in_file_path,count)
        count=count+1
        line=in_file.readline()
    out_file.close()
    in_file.close()
def split_gene_Promoter_file_from_dir(in_file,out_dir,sep,append=".dat"):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for chr_no in chr_no_list:
        split_gene_Promoter_file(chr_no,in_file,out_dir+"/chr"+chr_no+append,sep,append)
def split_gene_Promoter_file(chr_no,in_file_path,out_file_path,sep,append=".dat"):
    in_file=open(in_file_path,'r')
    out_file= open(out_file_path,'w')
    chr_no=str(chr_no)
    line=in_file.readline()
    count=0
    while line:
        count=count+1
        if count % 10000==0:
            print "chr%s promoter_split: %d of %s was processed" %(chr_no,count,in_file_path)
        match = re.search(r'(\d+|X|Y),(\d+),(\d+),([-|+]),([^\n]+)',line)
        chr_no_matched=match.group(1)
        if match and chr_no==chr_no_matched:
            promoter_S=match.group(2)
            promoter_T=match.group(3)
            strand=match.group(4)
            gene_name=match.group(5)
            out_file.write(chr_no_matched+sep+promoter_S+sep+promoter_T+sep+strand+sep+gene_name+"\n")
        line=in_file.readline()
    out_file.close()
    in_file.close()
# def statistic_cpg_count_from_dir(speices_name,seq_dir,promoter_splited_dir,in_append,out_dir,sep,append=".dat"):
#     if not os.path.exists(out_dir):
#         os.makedirs(out_dir)
#     for chr_no in chr_no_list:
#         statistic_cpg_count(chr_no,seq_dir+os.sep+speices_name+"dna.chromosome."+chr_no+".fa",promoter_splited_dir+os.sep+chr_no+in_append,out_dir+"/chr"+chr_no+append,sep,append)
    #根据序列文件以及输入的启动子范围文件,统计出每个启动子内的CpG位点数量.输出位置信息和CpG统计信息到out文件中
# def statistic_cpg_count(chr_no,seq_file_path,promoter_file_path,out_file_path,sep,append=".dat"):
#     seq_file = open(seq_file_path,'r')
#     promoter_file = open(promoter_file_path,'r')
#     out_file = open(out_file_path,'w')
#     line = promoter_file.readline()
#     match = re.search(r'([\d+|X|Y])'+sep+'(\d+)'+sep+'(\d+)'+sep+'([-|+])'+sep+'([^\s]+)',line)
#     promter_start=int(match.group(2))
#     promter_end=int(match.group(3))
#     strand=match.group(4)
#     gene_name=match.group(5)
#     out_str=""
#
#     start=time.clock()
#     line_seq=seq_file.readline()
#     match = re.search(r'chromosome:([^:]+):(\d*):(\d*):(\d*)',line_seq)
#     length=int(math.ceil(int(match.group(4))/60.0+1.0))
#     i=1
#     percent_100=length/100
#     now_pos=0
#     is_statisticing=False
#     promoter_cpg_now_count=0
#     is_end_with_c=False
#     while i < length+1:
#         line_seq = seq_file.readline()
#         if (i % percent_100== 0):
#             percent = i/percent_100
#             print "%d percent for statistic CpG of %s" % (percent,promoter_file.name)
#         #一行
#         if now_pos+60>promter_start and now_pos+60 < promter_end and not is_statisticing:
#             for j in range(promter_start-now_pos,60):
#                 if now_pos+j>=promter_start:
#                     if line[j-1]=='C' and line[j]=='G':
#                         promoter_cpg_now_count=promoter_cpg_now_count+1
#                     if j==59 and line[j]=='C':
#                         is_end_with_c
#             is_statisticing=True
#         elif now_pos+60>promter_start and now_pos+60 >= promter_end and not is_statisticing:
#
#         #把剩下的读完
#         elif now_pos+60>promter_end and is_statisticing:
#             is_statisticing=False
#         elif is_statisticing:
#         elif now_pos+60<promter_start:
#             pass
#         now_pos=now_pos+60

if __name__ == '__main__':
    sample_name="GSM1386021_2cell_mc_CG_paternal_plus"

    #将一个原始的bed文件切割成按照染色体划分的
    #batch_bed_to_methy(sample_name+".bed",sample_name)

    #准备文件夹和配置项名称
    data_dir="data"
    temp_dir_name=sample_name+os.sep+"tmp"
    stat_dir=sample_name+os.sep+"stat"
    out_dir=sample_name+os.sep+"ck_out"
    BASE_DIR = os.path.dirname(__file__)
    out_root = os.path.join(BASE_DIR, "out")

    species_name="Mus_musculus.GRCm38.dna.chromosome"

    if not os.path.exists(sample_name):
        os.makedirs(sample_name)
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    if not os.path.exists(temp_dir_name):
        os.makedirs(temp_dir_name)
    if not os.path.exists(stat_dir):
        os.makedirs(stat_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    #遍历染色体编号
    #for chr_no in chr_no_list:
    #!!!!!!!!!!!!!!!需要改这个的输出
    outer_out_dir="out"
    out_output=outer_out_dir+os.sep+"cg_count_40_plus"
    if not os.path.exists(outer_out_dir):
        os.makedirs(outer_out_dir)
    if not os.path.exists(out_output):
        os.makedirs(out_output)
    #!!!!!!!!!!!!!!!需要改成输入的文件
    interval_file_name="interval_40_plus.csv"
    out_interval_final_path=outer_out_dir+os.sep+interval_file_name
    filelist = os.listdir(out_root)
    for name in filelist:
        if not name.endswith(".csv"):
            path_to_remove=out_root +os.sep+ name
            try:
                os.remove(path_to_remove)  # delete files in path "directory"#
            except BaseException, e:
                print(str(e))
                print "this is a folder not a file!"
                #shutil.rmtree(path_to_remove) # delete an entire directory tree.
                print "files delete complete!"
    #infile_为筛选出的区间_cpg_区间,进行甲基化C_k曲线统计
    in_file = open(out_interval_final_path,'r')
    chr_no="1"
    line=in_file.readline()
    count_line=1
    seq=","
    window_length=200
    k_list=range(0,21)
    k_list2=[25,30,40,50,60,70,80,90,100,150,200,250,300]
    k_list.extend(k_list2)

    methy_count=0
    while line:
        pattern="(\d+)"+seq+"(\d+)"+seq+"(\d+)"+seq+"(\d+)\n"
        match = re.search(pattern,line)
        if match:
            start_site=int(match.group(3))
            next_site=int(match.group(4))
            end_site=next_site+window_length-1
            #扫描原始序列,并对比.bed文件,得出中间文件tmp/methylation_output_chr_no.txt,文件内容为"位置  碱基  甲基化水平"  统计文件为:stat/stat_chr_no.txt  里面内容包括:"总碱基数量,总甲基化水平,总甲基化碱基数,平均甲基化水平,所费时间"
            all_base_value,all_base_count,average_value=scan(data_dir+os.sep+species_name+"."+chr_no+".fa",sample_name+os.sep+"chr"+chr_no+".bed",temp_dir_name+os.sep+"methylation_output_"+chr_no+".txt",stat_dir+os.sep+"stat_"+chr_no+".txt",start_site,end_site,",")

            # #从文件中得出统计水平
            #in_file_stat= open(stat_dir+os.sep+"stat_"+chr_no+".txt",'r')
            # #总甲基化水平
            #all_base_value=float(in_file_stat.readline())

            if all_base_value!=0.0:

                methy_count=methy_count+1

            count_line=count_line+1

            print "line now:%d" % count_line
            # #总碱基数量
            # all_base_count=int(in_file.readline())
            # #平均甲基化水平
            # average_value=float(in_file.readline())


            #需要统计的k值list
            #k_list=[3,4,8]
            str_to_write_in="chr"+chr_no+"_"+str(start_site)+"_"+str(end_site)+""
            #遍历k值
            for k in k_list:
                #输出C(k)统计信息的文件
                out_file = open(out_output+os.sep+str_to_write_in+".txt",'a')
                #临时文件,即碱基和甲基化水平文件
                stat_file_path=temp_dir_name+os.sep+"methylation_output_"+chr_no+".txt"

                stat_file1=open(stat_file_path,'r')
                stat_file2=open(stat_file_path,'r')

                #统计k时的C(k),sum1=sum_{i=1}^{N-k}(a_i-\overline{a})(a_{i+k}-\overline{a}),sum2=sum_{i=1}^N(a_i-\overline{a})^2
                calculate(k,average_value,all_base_count,stat_file1,stat_file2,out_file,",",",")

                out_file.close()

                stat_file1.close()
                stat_file2.close()
                print "finish process %d: chr%s of %s" % (count_line,chr_no,sample_name)
        count_line=count_line+1
        os.remove(temp_dir_name+os.sep+"methylation_output_"+chr_no+".txt")
        #print "line %d" %count_line
        line=in_file.readline()
    in_file.close()

    print "finish ALL!%d" % methy_count
    #距离统计文件夹
    # distance_out_name="distance_out"
    # distribution_dir_name="count_distribute"
    # dis_info_dir_name="dis_info"
    # distinct_dis_info_dir_name="distinct_dis"
    #
    # dis_info_dir=sample_name+os.sep+dis_info_dir_name
    # distinct_dir=sample_name+os.sep+distinct_dis_info_dir_name
    # if not os.path.exists(dis_info_dir):
    #     os.makedirs(dis_info_dir)
    # if not os.path.exists(distinct_dir):
    #     os.makedirs(distinct_dir)
    #
    # chr_no='1'
    # append='.dat'
    # pattern=r'(\d+)\t([\d]+\.[\d]*)(\s)'
    # #输入.bed文件,根据输入的距离过滤条件来找出满足甲基化的前后两个甲基化的CpG位置和距离信息
    # calc_distance_filter_and_summary(sample_name+os.sep+"chr"+chr_no+".bed",dis_info_dir+os.sep+"chr"+chr_no+"_dis_info"+append,"<100",chr_no,pattern,",")
    #
    # distinct_dis_info(dis_info_dir+os.sep+"chr"+chr_no+"_dis_info"+append,distinct_dir+os.sep+"chr"+chr_no+"_distinct"+".csv")
    # extract_position_from_dir(sample_name,sample_name+os.sep+distance_out_name)
    # count_methylation_distribute_from_dir(sample_name+os.sep+distance_out_name,sample_name+os.sep+distribution_dir_name)
    # merge_count_distribute_from_dir(sample_name+os.sep+distribution_dir_name,sample_name+os.sep+distribution_dir_name+os.sep+"merged_count.csv")
    # get_all_gene_promoters("data"+os.sep+"Mus_musculus.GRCm38.83.chr.gtf",sample_name+os.sep+"Mus_musculus.GRCm38.83.chr.dat",1000,1000,",")
    # split_gene_Promoter_file_from_dir(sample_name+os.sep+"Mus_musculus.GRCm38.83.chr.dat",sample_name+os.sep+"promoter_splited",",",".dat")