import json
import itertools

configfile: "config.yaml"
R = config["R"]

# wildcards
DAT = glob_wildcards("code/00-get_dat-{x}.R").x
SIM = glob_wildcards("code/00-get_sim-{x}.R").x
SIZE = glob_wildcards("code/02-pbDEA-{x}.R").x

BIN = ["PAS","combined","asNA","sig"]
COR = ["gBH.LSL","binwise","IHW","gBH.TST"]
LOC = ["loc", "glb"]

STA = ["F1"]



# output
sim_dat = "data/sim/"; dat_dat = "data/dat/"
sim_out = "outs/sim/"; dat_out = "outs/dat/"

dat = expand(dat_dat+"00-raw/{dat}.rds", dat=DAT)
sim = expand(sim_dat+"00-raw/{sim}.rds", sim=SIM)
truth = expand(sim_dat+"00-truth/{sim}.rds", sim=SIM)
bulkDEA = expand(sim_out + "bulkDEA/{sim}.rds", sim=SIM)
pbDEA = expand(sim_out + "pbDEA/splatter,{size}.rds", sim=SIM, size=SIZE) + [sim_out + "pbDEA/muscat_LPS,12_2v2.rds"]

bbhw = expand(sim_out + "bbhw/splatter,{size},{bin},{cor},{loc}.rds",
              sim=SIM, size=SIZE, bin=BIN, cor=COR, loc=LOC) \
     + expand(sim_out + "bbhw/muscat_LPS,12_2v2,{bin},{cor},{loc}.rds",
              bin=BIN, cor=COR, loc=LOC)

sta = (
    expand(sim_out + "sta/{sta},splatter,{size},{bin},{cor},{loc}.rds",
           sta=STA, sim=["splatter"], size=SIZE,
           bin=BIN, cor=COR, loc=LOC)
    +
    expand(sim_out + "sta/{sta},muscat_LPS,12_2v2,{bin},{cor},{loc}.rds",
           sta=STA, sim=["muscat_LPS"],
           bin=BIN, cor=COR, loc=LOC)
)



# results
sim_res = {
    "sim": sim,
    "truth": truth,
    "bulkDEA": bulkDEA,
    "pbDEA": pbDEA,
    "bbhw": bbhw,
    "sta": sta
}

# visualization
VAL = sim_res.keys()

plt = []
for val in VAL:
    x = glob_wildcards("code/plt-"+val+"_{x}.R").x
    plt += expand("plts/{val}-{plt}.pdf", val=val, plt=x)

# SETUP ========================================================================
rule all: 
    input:
        "session_info.txt",
        [x for x in sim_res.values()], #plt,
        #[x for x in dat_res.values()], #qlt

        "session_info.txt"
        
rule session_info:
    priority: 100
    input:  "code/10-session_info.R"
    output: "session_info.txt"
    log:    "logs/session_info.Rout"
    shell:  '''
    {R} CMD BATCH --no-restore --no-save\
    "--args {output}" {input} {log}'''

# SIMULATION ===================================================================
rule get_sim:
    priority: 99
    input:  "code/00-get_sim.R",
            "code/00-get_sim-{sim}.R"
    output: 
        sim_dat+"00-raw/{sim}.rds",
        sim_dat+"00-truth/{sim}.rds"
    log:	"logs/get_dat-{sim}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        {input[1]} sce={output[0]} truth={output[1]}" {input[0]} {log}'''  

rule sim_bulkDEA:
    priority: 98
    input:  
        "code/01-bulkDEA.R",
        sim_dat+"00-raw/{sim}.rds"
    output: 
        sim_out + "bulkDEA/{sim}.rds"
    log:	"logs/sim_bulkDEA-{sim}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        {input[1]} {output[0]}" {input[0]} {log}'''  



rule sim_pbDEA:
    priority: 98
    input:  
        "code/02-pbDEA.R",
        "code/02-pbDEA-{size}.R",
        sim_dat+"00-raw/{sim}.rds"
    output: 
        sim_out + "pbDEA/{sim},{size}.rds"
    log:	"logs/sim_pbDEA-{sim},{size}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        {input[1]} {input[2]} {output[0]}" {input[0]} {log}'''  

rule sim_bbhw:
    priority: 97
    input:  
        "code/03-bbhw.R",
        sce=sim_dat+"00-raw/{sim}.rds",
        bulk=rules.sim_bulkDEA.output,
        pb=rules.sim_pbDEA.output
    output: 
        sim_out + "bbhw/{sim},{size},{bin},{cor},{loc}.rds"
    log:	"logs/sim_bbhw-{sim},{size},{bin},{cor},{loc}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        sce={input.sce} bulk={input.bulk}  pb={input.pb} res={output[0]}" {input[0]} {log}'''  


rule sim_sta:
    priority: 96
    input:
        "code/04-sta.R",
        "code/04-sta-{sta}.R",
        bbhw=sim_out + "bbhw/{sim},{size},{bin},{cor},{loc}.rds",
        truth=sim_dat+"00-truth/{sim}.rds"
    output:
        sim_out + "sta/{sta},{sim},{size},{bin},{cor},{loc}.rds"
    log:	"logs/sim_bbhw-{sta},{sim},{size},{bin},{cor},{loc}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        fun={input[1]} res={input.bbhw} truth={input.truth} sta={output[0]}" {input[0]} {log}'''  

############### Visualization #################
# for val in VAL:
#     rule:
#         priority: 90
#         input:  expand("code/plt-{val}_{{plt}}.R", val=val), x=res[val]
#         params: lambda wc, input: ";".join(input.x)
#         output: expand("plts/{val}-{{plt}}.pdf", val=val)
#         log:    expand("logs/plt_{val}-{{plt}}.Rout", val=val)
#         shell:  '''
#             {R} CMD BATCH --no-restore --no-save "--args\
#             {params} {output[0]}" {input[0]} {log}'''