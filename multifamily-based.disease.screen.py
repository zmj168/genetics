#!/usr/bin/python3
# -*- coding:utf-8 -*-
####################################################################
#
####################################################################
'''
\033[1;34;47mUSAGE:
      Description: 基于多家系的复杂疾病位点筛选，本程序使用CMsiteOnly_Case (共有特有分析), Exomiser (基于临床表型的致病性分析), pVAAST(连锁和相关性分析)和FARVAT (family-based gene-burden分析)，并使用Intervar进行了致病性预测，整合上述方法，筛选致病性基因和位点，并对其进行DisGenet数据库注释
      分析模块:1 InterVar, 2 CMsiteOnly_Case, 3 Exomiser, 4 pVAAST, 5, FARVAT, 6 Candigene_integration      
      模块间依赖关系:  relation = {
                                    InterVar:[]
                                    CMsiteOnly_Case:[InterVar]
                                    Exomiser:[]
                                    pVAAST:[]
                                    FARVAT:[]
                                    Candigene_integration:[Intervar,CMsiteOnly_Case,Exomiser,pVAAST,FARVAT]
                                   }                    
      To run : python3 multifamily-based.disease.screen.py [options]
      options: -c config file [required]
      e.g.:
            python3 multifamily-based.disease.screen.py -c multifamily-based.disease.screen.config.example.ini\033[0m
'''

Version = "v1.0.0"
print ("\033[1;34;47mmultifamily-based.disease.screen pipeline{}Version: {}.\033[0m".format("\n",Version))

import os
import sys
import re
import glob

import argparse
from argparse import RawTextHelpFormatter

Bin=os.path.abspath(os.path.dirname(__file__))
sys.path.append(Bin+'/../lib')
from PipMethod import myconf,generateShell,mkdir

__author__='zhangshouwei'
__mail__='shouweizhang@genome.cn'

def check_file(file):
    if not os.path.isfile(file):
        sys.stderr.write('\033[1;31;40mmultifamily-based.disease.screen pipeline - ERROR - input: No such file {0}. Please ensure the correct file path'.format(file))
        exit(1)

def check_path(path,status):
    """status: new: mkdir
               quit: exit
    """
    if not os.path.exists(path):
        if status == 'quit':
            sys.stderr.write('\033[1;31;40mmultifamily-based.disease.screen pipeline - ERROR - input: No such path {0}. Please ensure the correct path'.format(path))
            exit(1)
        elif status == 'new':
            os.makedirs(path)

def generateShell_exe(shell,cmd,module):
    with open(shell,'w') as shell_open:
        shell_open.writelines('{}\n'.format(cmd))
    try:
        os.system(cmd)
    except:
        print ('\033[1;31;40mPlease check this {}!Some errors has occured!'.format(module))
        sys.exit(1)
    else:
        print ('\033[1;34;47m{} proceed successfully\033[0m'.format(module))
        print ('{}\n'.format(cmd))

def main():
    '''参数设定'''
    parser = argparse.ArgumentParser(description=__doc__,formatter_class=RawTextHelpFormatter,epilog='author:\t{0}\nmail:\t{1}\n'.format(__author__,__mail__))
    parser.add_argument('-c','--config',help="multifamily-based.disease.screen pipeline configuration file: eg: multifamily-based.disease.screen.config.example.ini.",dest='config',type=str,required=True)
    argv = vars(parser.parse_args())

    config = argv['config'].strip()
    allmodule = {'InterVar':0,
                'CMsiteOnly_Case':0,
                'Exomiser':0,
                'FARVAT':0,
                'pVAAST':0,
                'Candigene_integration':0}

    dependent = {'InterVar':[],
                 'CMsiteOnly_Case':['InterVar'],
                 'Exomiser':[],
                 'pVAAST':[],
                 'FARVAT':[],
                 'Candigene_integration':['InterVar','CMsiteOnly_Case','Exomiser','pVAAST','FARVAT']}

    cfg = myconf()
    cfg.readfp(open(config))
   
    analysis_content = cfg.get('sys','analysisModules').split(',')
    ## 检查填写的分析内容是否正确 ##
    for i in analysis_content:
        s = set(dependent[i])
        if not s.issubset(set(analysis_content)):
            print ('Error:\tyou need do %s before %s' %(', '.join(dependent[i]), i))
            exit(1)
        else:pass 

    ## 进一步检查分析模块,并添加分析标志 ##
    moduleNum = 0
    for i in analysis_content:
        if i in allmodule:
            allmodule[i] += 1
            moduleNum += 1
        else:
            print ('\033[1;31;40mWarning: the analysis content {} is wrong'.format(i))

    if moduleNum == 0:
        print ('\033[1;31;40mThere is no Analysis Modules in {}'.format(config))
        exit(1)

    odir = cfg.get('sys','odir')
    moduledir = cfg.get('sys','moduledir')
    monitorOptions = cfg.get('sys','monitorOptions')
    vcf2peddataset = cfg.get('sys','vcf2peddataset')
    familyhpo = cfg.get('sys','familyhpo')
    gedition = cfg.get('sys','gedition')

    check_file(vcf2peddataset);check_file(familyhpo)
    check_path(odir,'exit');check_path(moduledir,'exit')

    tasklog = os.path.join(odir,'tasklog')
    mkdir([tasklog])
    shelldir = os.path.join(odir,'shell')
    processdir = os.path.join(odir,'process')
    listdir = os.path.join(odir,'list')
    configdir = os.path.join(odir,'config')
    mkdir([shelldir,processdir,listdir,configdir])
    dependent_list = []

    ## Analysis Module Part ##
    ## Intervar Module ##
    if allmodule['InterVar'] == 1:
        intervar_config_find = glob.glob('{0}/InterVar/bin/InterVar.ini'.format(moduledir))
        if len(intervar_config_find) == 1:
            intervar_config = intervar_config_find[0]
        else:
            print ('\033[1;31;40mPlease chech whether InterVar config existed or not!')
            exit(1)
        
        ## config modify ##
        InterVar_cfg = myconf()
        InterVar_cfg.readfp(open(intervar_config))

        InterVar_cfg.set('software','annovar',cfg.get('software','annovar'))
        InterVar_cfg.set('software','python3',cfg.get('software','python3'))
        InterVar_cfg.set('software','intervar',cfg.get('software','intervar'))
        InterVar_cfg.set('software','bgzip',cfg.get('software','bgzip'))
        InterVar_cfg.set('software','tabix',cfg.get('software','tabix'))
        InterVar_cfg.set('software','monitor',cfg.get('software','monitor'))
        InterVar_cfg.set('database','database_intervar',cfg.get('database','database_intervar'))
        InterVar_cfg.set('database','database_locat',cfg.get('database','database_locat'))
        InterVar_cfg.set('memory','Intervar_run',cfg.get('InterVar_memory','Intervar_run'))
        InterVar_cfg.set('memory','family_integrate',cfg.get('InterVar_memory','family_integrate'))
        InterVar_cfg.set('memory','add_genotype',cfg.get('InterVar_memory','add_genotype'))
        InterVar_cfg.set('memory','all_intergrate',cfg.get('InterVar_memory','all_intergrate'))
        InterVar_cfg.set('memory','upload',cfg.get('InterVar_memory','upload'))
        InterVar_cfg.set('genomeedition','gedition',cfg.get('sys','gedition'))
        InterVar_cfg.write(open('{0}/InterVar.ini'.format(configdir), "w"))

        shell = os.path.join(tasklog,'InterVar.sh')
        cmd = '{0} {1}/InterVar/bin/InterVar_pipeline.py -ds {2} -conf {3}/InterVar.ini -moption "taskmonitor -q {4}" -md InterVar -o {5}'.format(cfg.get('software','python3'),moduledir,vcf2peddataset,configdir,cfg.get('sys','queue'),odir)
        generateShell_exe(shell,cmd,'InterVar')
        dependent_list.append('%s/InterVar_dependence.txt' %(listdir))
    else:
        print ("\033[1;34;47mPlease check pipeline config carefully! You must proceed {} module\033[0m".format('InterVar'))
        exit(1)

    ## CMsiteOnly_Case Module ##
    if allmodule['CMsiteOnly_Case'] == 1:
        cmsc_config_find = glob.glob('{0}/CMsiteOnly_Case/bin/CMsiteOnly_Case.ini'.format(moduledir))
        if len(cmsc_config_find) == 1:
            cmsc_config = cmsc_config_find[0]
        else:
            print ('\033[1;31;40mPlease chech whether CMsiteOnly_Case config existed or not!')
            exit(1)
        
        ## config modify ##
        cmsc_cfg = myconf()
        cmsc_cfg.readfp(open(cmsc_config))

        cmsc_cfg.set('software','python3',cfg.get('software','python3'))
        cmsc_cfg.set('software','perl',cfg.get('software','perl'))
        cmsc_cfg.set('software','monitor',cfg.get('software','monitor'))
        cmsc_cfg.set('memory','rarebasic',cfg.get('CMsiteOnly_Case_memory','rarebasic'))
        cmsc_cfg.set('memory','familyscreen',cfg.get('CMsiteOnly_Case_memory','familyscreen'))
        cmsc_cfg.set('memory','pathogenityscreen',cfg.get('CMsiteOnly_Case_memory','pathogenityscreen'))
        cmsc_cfg.set('memory','gnomADscreen',cfg.get('CMsiteOnly_Case_memory','gnomADscreen'))
        cmsc_cfg.set('memory','stat',cfg.get('CMsiteOnly_Case_memory','stat'))
        cmsc_cfg.set('memory','upload',cfg.get('CMsiteOnly_Case_memory','upload'))
        cmsc_cfg.set('genomeedition','gedition',cfg.get('sys','gedition')) 
        cmsc_cfg.set('para','1000g',cfg.get('CMsiteOnly_Case_para','1000g'))
        cmsc_cfg.set('para','gnomAD',cfg.get('CMsiteOnly_Case_para','gnomAD'))
        cmsc_cfg.set('para','minfamily',cfg.get('CMsiteOnly_Case_para','minfamily'))
        cmsc_cfg.write(open('{0}/CMsiteOnly_Case.ini'.format(configdir), "w"))

        shell = os.path.join(tasklog,'CMsiteOnly_Case.sh')
        cmd = '{0} {1}/CMsiteOnly_Case/bin/CMsiteOnly_Case.pipeline.py -ds {2} -conf {3}/CMsiteOnly_Case.ini -il {4}/{5}.list -moption "taskmonitor -q {6}" -md CMsiteOnly_Case -o {7}'.format(cfg.get('software','python3'),moduledir,vcf2peddataset,configdir,listdir,'InterVar',cfg.get('sys','queue'),odir)
        generateShell_exe(shell,cmd,'CMsiteOnly_Case')
        dependent_list.append('%s/CMsiteOnly_Case_dependence.txt' %(listdir))
    else:
        print ("\033[1;34;47mPlease check pipeline config carefully! You must proceed {} module\033[0m".format('CMsiteOnly_Case'))
        exit(1)

    ## Exomiser Module ##
    if allmodule['Exomiser'] == 1:
        Exomiser_config_find = glob.glob('{0}/Exomiser/bin/Exomiser.ini'.format(moduledir))
        if len(Exomiser_config_find) == 1:
            Exomiser_config = Exomiser_config_find[0]
        else:
            print ('\033[1;31;40mPlease chech whether Exomiser config existed or not!')
            exit(1)

        ## config modify ##
        Exomiser_cfg = myconf()
        Exomiser_cfg.readfp(open(Exomiser_config))

        Exomiser_cfg.set('software','python3',cfg.get('software','python3'))
        Exomiser_cfg.set('software','java',cfg.get('software','java'))
        Exomiser_cfg.set('software','monitor',cfg.get('software','monitor'))
        Exomiser_cfg.set('software','exomiser',cfg.get('software','exomiser'))
        Exomiser_cfg.set('memory','exomconfig',cfg.get('Exomiser_memory','exomconfig'))
        Exomiser_cfg.set('memory','exomrun',cfg.get('Exomiser_memory','exomrun'))
        Exomiser_cfg.set('memory','exomfilter',cfg.get('Exomiser_memory','exomfilter'))
        Exomiser_cfg.set('memory','stat',cfg.get('Exomiser_memory','stat'))
        Exomiser_cfg.set('memory','upload',cfg.get('Exomiser_memory','upload'))
        Exomiser_cfg.set('genomeedition','gedition',cfg.get('sys','gedition'))
        Exomiser_cfg.set('database','exomiser_data',cfg.get('database','exomiser_data'))
        Exomiser_cfg.set('database','version',cfg.get('database','version'))
        Exomiser_cfg.set('para','minfamily',cfg.get('Exomiser_para','minfamily'))
        Exomiser_cfg.write(open('{0}/Exomiser.ini'.format(configdir), "w"))

        shell = os.path.join(tasklog,'Exomiser.sh')
        cmd = '{0} {1}/Exomiser/bin/Exomiser.pipeline.py -ds {2} -conf {3}/Exomiser.ini -fhpo {4} -ih {5} -moption "taskmonitor -q {6}" -md Exomiser -o {7}'.format(cfg.get('software','python3'),moduledir,vcf2peddataset,configdir,cfg.get('sys','familyhpo'),cfg.get('Exomiser_para','inheritance'),cfg.get('sys','queue'),odir)
        generateShell_exe(shell,cmd,'Exomiser')
        dependent_list.append('%s/Exomiser_dependence.txt' %(listdir))
    else:
        print ("\033[1;34;47mPlease check pipeline config carefully! You must proceed {} module\033[0m".format('Exomiser'))
        exit(1)

    ## FARVAT Module ##
    if allmodule['FARVAT'] == 1:
        FARVAT_config_find = glob.glob('{0}/FARVAT/bin/FARVAT.ini'.format(moduledir))
        if len(FARVAT_config_find) == 1:
            FARVAT_config = FARVAT_config_find[0]
        else:
            print ('\033[1;31;40mPlease chech whether FARVAT config existed or not!')
            exit(1)

        ## config modify ##
        FARVAT_cfg = myconf()
        FARVAT_cfg.readfp(open(FARVAT_config))
        FARVAT_cfg.set('software','python3',cfg.get('software','python3'))
        FARVAT_cfg.set('software','monitor',cfg.get('software','monitor'))
        FARVAT_cfg.set('software','bgzip',cfg.get('software','bgzip'))
        FARVAT_cfg.set('software','tabix',cfg.get('software','tabix'))
        FARVAT_cfg.set('software','vcftools',cfg.get('software','vcftools'))
        FARVAT_cfg.set('software','vcfmerge',cfg.get('software','vcfmerge'))
        FARVAT_cfg.set('memory','prepedvcf',cfg.get('FARVAT_memory','prepedvcf'))
        FARVAT_cfg.set('memory','snpvcfmerge',cfg.get('FARVAT_memory','snpvcfmerge'))
        FARVAT_cfg.set('memory','farvat',cfg.get('FARVAT_memory','farvat'))
        FARVAT_cfg.set('memory','stat',cfg.get('FARVAT_memory','stat'))
        FARVAT_cfg.set('memory','upload',cfg.get('FARVAT_memory','upload'))
        FARVAT_cfg.set('genomeedition','gedition',cfg.get('sys','gedition'))
        FARVAT_cfg.set('database','1000g',cfg.get('database','1000g'))
        FARVAT_cfg.set('database','geneinfo',cfg.get('database','geneinfo'))
        FARVAT_cfg.set('para','rare',cfg.get('FARVAT_para','rare'))
        FARVAT_cfg.set('para','common',cfg.get('FARVAT_para','common'))
        FARVAT_cfg.set('para','pvalue',cfg.get('FARVAT_para','pvalue'))
        FARVAT_cfg.write(open('{0}/FARVAT.ini'.format(configdir), "w"))

        shell = os.path.join(tasklog,'FARVAT.sh')
        cmd = '{0} {1}/FARVAT/bin/FARVAT.pipeline.py -ds {2} -conf {3}/FARVAT.ini -moption "taskmonitor -q {4}" -md FARVAT -o {5}'.format(cfg.get('software','python3'),moduledir,vcf2peddataset,configdir,cfg.get('sys','queue'),odir)
        generateShell_exe(shell,cmd,'FARVAT')
        dependent_list.append('%s/FARVAT_dependence.txt' %(listdir))
    else:
        print ("\033[1;34;47mPlease check pipeline config carefully! You must proceed {} module\033[0m".format('FARVAT'))
        exit(1)

    ## pVAAST Module ##
    if allmodule['pVAAST'] == 1:
        pVAAST_config_find = glob.glob('{0}/pVAAST/bin/pVAAST.ini'.format(moduledir))
        if len(pVAAST_config_find) == 1:
            pVAAST_config = pVAAST_config_find[0]
        else:
            print ('\033[1;31;40mPlease chech whether pVAAST config existed or not!')
            exit(1)

        ## config modify ##
        pVAAST_cfg = myconf()
        pVAAST_cfg.readfp(open(pVAAST_config))
        pVAAST_cfg.set('software','python3',cfg.get('software','python3'))
        pVAAST_cfg.set('software','monitor',cfg.get('software','monitor'))
        pVAAST_cfg.set('software','perl',cfg.get('software','perl'))
        pVAAST_cfg.set('software','vcf2cdr',cfg.get('software','vcf2cdr'))
        pVAAST_cfg.set('software','VAAST',cfg.get('software','VAAST'))
        pVAAST_cfg.set('memory','vcftocdr',cfg.get('pVAAST_memory','vcftocdr'))
        pVAAST_cfg.set('memory','pVAASTConfig',cfg.get('pVAAST_memory','pVAASTConfig'))
        pVAAST_cfg.set('memory','pVAAST1',cfg.get('pVAAST_memory','pVAAST1'))
        pVAAST_cfg.set('memory','genefocus',cfg.get('pVAAST_memory','genefocus'))
        pVAAST_cfg.set('memory','pVAAST2',cfg.get('pVAAST_memory','pVAAST2'))
        pVAAST_cfg.set('memory','stat',cfg.get('pVAAST_memory','stat'))
        pVAAST_cfg.set('memory','upload',cfg.get('pVAAST_memory','upload'))
        pVAAST_cfg.set('genomeedition','gedition',cfg.get('sys','gedition'))
        pVAAST_cfg.set('database','genomefa',cfg.get('database','genomefa'))
        pVAAST_cfg.set('database','gff3',cfg.get('database','gff3'))
        pVAAST_cfg.set('para','pVAASTgw',cfg.get('pVAAST_para','pVAASTgw'))
        pVAAST_cfg.set('para','pVAASTd',cfg.get('pVAAST_para','pVAASTd'))
        pVAAST_cfg.set('para','lod',cfg.get('pVAAST_para','lod'))
        pVAAST_cfg.set('para','pvalue',cfg.get('pVAAST_para','pvalue'))
        pVAAST_cfg.write(open('{0}/pVAAST.ini'.format(configdir), "w"))

        shell = os.path.join(tasklog,'pVAAST.sh')
        cmd = '{0} {1}/pVAAST/bin/pVAAST.pipeline.py -ds {2} -conf {3}/pVAAST.ini -moption "taskmonitor -q {4}" -ih {5} -md pVAAST -o {6}'.format(cfg.get('software','python3'),moduledir,vcf2peddataset,configdir,cfg.get('sys','queue'),cfg.get('pVAAST_para','inheritance'),odir)
        generateShell_exe(shell,cmd,'pVAAST')
        dependent_list.append('%s/pVAAST_dependence.txt' %(listdir))
    else:
        print ("\033[1;34;47mPlease check pipeline config carefully! You must proceed {} module\033[0m".format('pVAAST'))
        exit(1)

    ## Candigene_integration Module ##
    if allmodule['Candigene_integration'] == 1:
        cdi_config_find = glob.glob('{0}/Candigene_integration/bin/Candigene_integration.ini'.format(moduledir))
        if len(cdi_config_find) == 1:
            cdi_config = cdi_config_find[0]
        else:
            print ('\033[1;31;40mPlease chech whether Candigene_integration config existed or not!')
            exit(1)

        ## config modify ##
        cdi_cfg = myconf()
        cdi_cfg.readfp(open(cdi_config))
        cdi_cfg.set('software','python3',cfg.get('software','python3'))
        cdi_cfg.set('software','monitor',cfg.get('software','monitor'))
        cdi_cfg.set('memory','candigenescreen',cfg.get('Candigene_integration_memory','candigenescreen'))
        cdi_cfg.set('memory','gene2entrezid-intervar',cfg.get('Candigene_integration_memory','gene2entrezid-intervar'))
        cdi_cfg.set('memory','digGenetanno',cfg.get('Candigene_integration_memory','digGenetanno'))
        cdi_cfg.set('memory','upload',cfg.get('Candigene_integration_memory','upload'))
        cdi_cfg.set('genomeedition','gedition',cfg.get('sys','gedition'))
        cdi_cfg.set('database','diGenent',cfg.get('database','diGenent'))
        cdi_cfg.set('database','gene2ensembl',cfg.get('database','gene2ensembl'))
        cdi_cfg.set('para','minmethod',cfg.get('Candigene_integration_para','minmethod'))
        cdi_cfg.set('para','1000g',cfg.get('Candigene_integration_para','1000g'))
        cdi_cfg.set('para','gnomAD',cfg.get('Candigene_integration_para','gnomAD'))
        cdi_cfg.write(open('{0}/Candigene_integration.ini'.format(configdir), "w"))
 
        ## Methods results list prepare ##
        intervar_list = '{0}/{1}.final.list'.format(listdir,'InterVar')
        cmsc_list = '{0}/{1}.final.list'.format(listdir,'CMsiteOnly_Case')
        Exomiser_list = '{0}/{1}.final.list'.format(listdir,'Exomiser')
        FARVAT_list = '{0}/{1}.final.list'.format(listdir,'FARVAT')
        pVAAST_list = '{0}/{1}.final.list'.format(listdir,'pVAAST')
        with open('{0}/allmethods.list'.format(listdir),'w') as ON:
            ON.writelines('CMsiteOnly_Case\t%s\n' %(cmsc_list))
            ON.writelines('Exomiser\t%s\n' %(Exomiser_list))
            ON.writelines('FARVAT\t%s\n' %(FARVAT_list))
            ON.writelines('pVAAST\t%s\n' %(pVAAST_list))

        shell = os.path.join(tasklog,'Candigene_integration.sh')
        cmd = '{0} {1}/Candigene_integration/bin/Candigene_integration.pipeline.py -mm {2}/allmethods.list -iv {3} -conf {4}/Candigene_integration.ini -moption "taskmonitor -q {5}" -md Candigene_integration -o {6}'.format(cfg.get('software','python3'),moduledir,listdir,intervar_list,configdir,cfg.get('sys','queue'),odir)
        generateShell_exe(shell,cmd,'Candigene_integration')
        dependent_list.append('%s/Candigene_integration_dependence.txt' %(listdir))
    else:
        print ("\033[1;34;47mPlease check pipeline config carefully! You must proceed {} module\033[0m".format('Candigene_integration'))
        exit(1)

    ## qsub shell ##
    for i in dependent_list:
        os.system('cat {0} >> {1}/all_dependence.txt'.format(i,listdir))

    with open('{}/multifamily-based.disease.screen_qsub.sh'.format(odir),'w') as qsub_sh:
        qsub_sh.writelines('{0} taskmonitor {1} -i {2}/all_dependence.txt'.format(cfg.get('software','monitor'),cfg.get('sys','monitorOptions'),listdir,listdir))
    print ('\033[1;34;47m\nIf you have any questions, please contact shouweizhang@genome.cn! Have a happy cooperation!\n')
    os.system('cd {0} && sh multifamily-based.disease.screen_qsub.sh'.format(odir))

if __name__ == '__main__':
    main()
