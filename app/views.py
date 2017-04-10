from app.models import CellLine, Chromosome, Gene, Exon, Transcript, Fusion
import json, csv
from neomodel import db

# Create your views here.
from django.http import HttpResponse

def search_for_cell_line(request,c_line):
    response = {}
    header = get_header()
    #get_cell_line_from_disease("Colon adenocarcinoma")    
    # recupero fusioni nella linea cellulare
    fusions = []
    #se per tutte le linee cellulari mostra tutte le fusioni
    if c_line == "ALL":
        for c_l in CellLine.nodes.all():
            for fusion in c_l.happen:
                fusions.append(fusion)
                #print(fusion)
    else:
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            fusions.append(fusion)
    #print(fusions)
    rows = build_rows(fusions)
    #print(rows)
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_chromosome(request,chromos1,chromos2,c_line):
    response = {}
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    #TUTTE LE LINEE CELLULARI, UN CROMOSOMA 
    if c_line == "ALL" and chromos1 != "ALL" and chromos2 == "ALL" :
        c = Chromosome.nodes.get(chromosome = chromos1)
        for fusion in c.fromFusiontoChromosome:
            fusions.append(fusion)
    #TUTTE LE LINEE CELLULARI, ENTRAMBI I CROMOSOMI
    elif c_line == "ALL" and chromos1 != "ALL" and chromos2 != "ALL":
        c1 = Chromosome.nodes.get(chromosome = chromos1)
        c2 = Chromosome.nodes.get(chromosome = chromos2)
        for fusion in c1.fromFusiontoChromosome:
            if c2 in fusion.at_chromosome:
                fusions.append(fusion)
        for fusion in c2.fromFusiontoChromosome:
            if c1 in fusion.at_chromosome:
                fusions.append(fusion)
    #LINEA CELLULARE SPECIFICATA, UN CROMOSOMA
    elif c_line != "ALL" and chromos1 != "ALL" and chromos2 == "ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion.at_chromosome.filter(chromosome__exact=chromos1):
                fusions.append(fusion)
    #LINEA CELLULARE SPECIFICATA, ENTRAMBI I CROMOSOMI SPECIFICATI
    elif c_line != "ALL" and chromos1 != "ALL" and chromos2 != "ALL":
        c1 = Chromosome.nodes.get(chromosome = chromos1)
        c2 = Chromosome.nodes.get(chromosome = chromos2)
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if c1 in fusion.at_chromosome and c2 in fusion.at_chromosome:
                fusions.append(fusion)
    
                    
    rows = build_rows(fusions)
    #print(rows)
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))
        
def OLD_search_for_chromosome(request,chromos,c_line):
    response = {}
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    #se ho specificato solo il cromosoma, cerco tutte le fusioni in tutte le linee cellulari che coinvolgono il cromosoma
    if c_line == "ALL" and chromos != "ALL":
        c = Chromosome.nodes.get(chromosome = chromos)
        for fusion in c.fromFusiontoChromosome:
            fusions.append(fusion)
    #se ho specificato solo la linea cellulare, cerco tutte le fusioni che avvengono in uqella linea cellulare nel determinato intervallo
    elif c_line != "ALL" and chromos == "ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            fusions.append(fusion)
    #se ho sia linea cellulare che cromosoma specificati, cerco tutte le fusioni nella linea cellulare che coinvolgono il cromosoma
    elif c_line != "ALL" and chromos != "ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion.at_chromosome.filter(chromosome__exact=chromos):
                fusions.append(fusion)
    
    rows = build_rows(fusions)
    #print(rows)
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_gene(request,gene_one,gene_two,c_line):
    response = {}
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    
    #se ho specificato il gene ma non la linea cellulare, mostro tutte le fusioni che coinvolgono il gene
    if c_line == "ALL" and gene_one != "ALL" and gene_two == "ALL":
        if "ENSG" in gene_one:
            g = Gene.nodes.get(ensid = gene_one)
            if g.had:
                for fusion in g.had:
                    fusions.append(fusion)
            if g.fromFusionToGene:
                for fusion in g.fromFusionToGene:
                    fusions.append(fusion)
        else:
            g = Gene.nodes.get(symbol = gene_one)
            if g.had:
                for fusion in g.had:
                    fusions.append(fusion)
            if g.fromFusionToGene:
                for fusion in g.fromFusionToGene:
                    fusions.append(fusion)
    elif c_line == "ALL" and gene_one != "ALL" and gene_two != "ALL":
        if "ENSG" in gene_one and "ENSG" in gene_two:
            g = Gene.nodes.get(ensid = gene_one)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].ensid == gene_two:
                        fusions.append(fusion)
        
            g = Gene.nodes.get(ensid = gene_two)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].ensid == gene_one:
                        fusions.append(fusion)
        elif "ENSG" in gene_one and "ENSG" not in gene_two:
            g = Gene.nodes.get(ensid = gene_one)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].symbol == gene_two:
                        fusions.append(fusion)
        
            g = Gene.nodes.get(symbol = gene_two)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].ensid == gene_one:
                        fusions.append(fusion)
        elif "ENSG" not in gene_one and "ENSG" in gene_two:
            g = Gene.nodes.get(symbol = gene_one)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].ensid == gene_two:
                        fusions.append(fusion)
        
            g = Gene.nodes.get(ensid = gene_two)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].symbol == gene_one:
                        fusions.append(fusion)
        else:
            g = Gene.nodes.get(symbol = gene_one)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].symbol == gene_two:
                        fusions.append(fusion)
                        
            g = Gene.nodes.get(symbol = gene_two)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].symbol == gene_one:
                        fusions.append(fusion)    
    #se ho specificato la linea cellulare ma non il gene, mostro tutte le fusioni per quella linea cellulare (ANALOGO A SEARCH FOR CELL_LINE)
    elif c_line != "ALL" and gene_one == "ALL" and gene_two == "ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            fusions.append(fusion)     
    #se ho specificato sia la linea cellulare che il gene, cerco tutte le fusioni che coinvolgono quel gene nella linea cellulare
    elif c_line != "ALL" and gene_one != "ALL" and gene_two == "ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion.fromGeneToFusion.filter(symbol__exact=gene_one) or fusion.fromGeneToFusion.filter(ensid__exact=gene_one) or fusion.with_gene.filter(symbol__exact=gene_one) or fusion.with_gene.filter(ensid__exact=gene_one): #PROBLEMA
                fusions.append(fusion)
    elif c_line != "ALL" and gene_one != "ALL" and gene_two != "ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if (fusion.fromGeneToFusion.filter(symbol__exact=gene_one) and fusion.with_gene.filter(symbol__exact=gene_two)) or (fusion.fromGeneToFusion.filter(symbol__exact=gene_two) and fusion.with_gene.filter(symbol__exact=gene_one)) or (fusion.fromGeneToFusion.filter(ensid__exact=gene_one) and fusion.with_gene.filter(ensid__exact=gene_two)) or (fusion.fromGeneToFusion.filter(ensid__exact=gene_two) and fusion.with_gene.filter(ensid__exact=gene_one)): #PROBLEMA
                fusions.append(fusion)
            
    rows = build_rows(fusions)
    #print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))



def search_for_exon(request,exon_one,exon_two,c_line):
    response = {}
    header = get_fc_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    
    #cell_line all, primo esone si, secondo esone no, tutte le fusioni che coinvolgono questo esone
    if c_line == "ALL" and exon_one!="ALL" and exon_two == "ALL":
        e = Exon.nodes.get(exon = exon_one)
        for fcfusion in e.fromFusionToExon:
            for fusion in fcfusion.fromFusionToFusionCatcher:
                fusions.append(fusion)
    #cell_line all, primo esone si, secondo esone si (011), tutte le fusioni che convolgono la coppia di esoni
    elif c_line == "ALL" and exon_one!="ALL" and exon_two!="ALL":
        e = Exon.nodes.get(exon = exon_one)
        for fcfusion in e.fromFusionToExon:
            if fcfusion.at_exon.filter(exon__exact=exon_two):
                fusions.append(fcfusion.fromFusionToFusionCatcher.all()[0])
    #cell_line si, primo esone no, secondo esone no, ANALOGO A SEARCH FOR CELL_LINE
    elif c_line!="ALL" and exon_one=="ALL" and exon_two=="ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            fusions.append(fusion)     
    #cell_line si, primo esone si, secondo esone no 
    elif c_line != "ALL" and exon_one != "ALL" and exon_two == "ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            for fcfusion in fusion.with_fc_script:
                if fcfusion.at_exon.filter(exon__exact=exon_one):
                    fusions.append(fusion)
    #cell_line si, primo esone si, secondo esone si 
    elif c_line != "ALL" and exon_one!="ALL" and exon_two != "ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            for fcfusion in fusion.with_fc_script:
                if fcfusion.at_exon.filter(exon__exact=exon_one) and fcfusion.at_exon.filter(exon__exact=exon_two):
                    fusions.append(fusion)

    rows = build_fc_rows(fusions)
    #print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_transcript(request,transcript_one,transcript_two,c_line):
    response = {}
    header = get_fc_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    
    #linea cellulare non specificata, un trascritto specificato
    if c_line=="ALL" and transcript_one!="ALL" and transcript_two=="ALL":
        for couple in Transcript.nodes.get(transcript = transcript_one).in_couple:
            for fcfusion in couple.fromFusionToCouple:
                fusions.append(fcfusion.fromFusionToFusionCatcher.all()[0])
        for couple in Transcript.nodes.get(transcript = transcript_one).fromCoupleToTranscript:
            for fcfusion in couple.fromFusionToCouple:
                fusions.append(fcfusion.fromFusionToFusionCatcher.all()[0])
    #linea cellulare non specificata, coppia di trascritti specificata
    elif c_line=="ALL" and transcript_one!="ALL" and transcript_two!="ALL":
        for couple in Transcript.nodes.get(transcript = transcript_one).in_couple:
            if couple.with_other_transcript.filter(transcript__exact=transcript_two) or couple.fromTranscriptToCouple.filter(transcript__exact=transcript_two):
                for fcfusion in couple.fromFusionToCouple:
                    fusions.append(fcfusion.fromFusionToFusionCatcher.all()[0])
        for couple in Transcript.nodes.get(transcript = transcript_one).fromCoupleToTranscript:
            if couple.with_other_transcript.filter(transcript__exact=transcript_two) or couple.fromTranscriptToCouple.filter(transcript__exact=transcript_two):
                for fcfusion in couple.fromFusionToCouple:
                    fusions.append(fcfusion.fromFusionToFusionCatcher.all()[0])            
        for couple in Transcript.nodes.get(transcript = transcript_two).in_couple:
            if couple.with_other_transcript.filter(transcript__exact=transcript_two) or couple.fromTranscriptToCouple.filter(transcript__exact=transcript_two):
                for fcfusion in couple.fromFusionToCouple:
                    fusions.append(fcfusion.fromFusionToFusionCatcher.all()[0])            
        for couple in Transcript.nodes.get(transcript = transcript_two).fromCoupleToTranscript:
            if couple.with_other_transcript.filter(transcript__exact=transcript_one) or couple.fromTranscriptToCouple.filter(transcript__exact=transcript_one):
                for fcfusion in couple.fromFusionToCouple:
                    fusions.append(fcfusion.fromFusionToFusionCatcher.all()[0])
    #linea cellulare specificata, un trascritto specificato:
    elif c_line!="ALL" and transcript_one!="ALL" and transcript_two=="ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            for fcfusion in fusion.with_fc_script:
                for couple in fcfusion.with_trans_couple:
                    if couple.fromTranscriptToCouple.filter(transcript__exact=transcript_one) or couple.with_other_transcript.filter(transcript__exact=transcript_one):
                        fusions.append(fusion)
    #linea cellulare specificata, coppia di trascritti specificata
    elif c_line!="ALL" and transcript_one!="ALL" and transcript_two!="ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            for fcfusion in fusion.with_fc_script:
                for couple in fcfusion.with_trans_couple:
                    if (couple.fromTranscriptToCouple.filter(transcript__exact=transcript_one) and couple.with_other_transcript.filter(transcript__exact=transcript_two)) or (couple.fromTranscriptToCouple.filter(transcript__exact=transcript_two) and couple.with_other_transcript.filter(transcript__exact=transcript_one)):
                        fusions.append(fusion)

    rows = build_fc_rows(fusions)
    #print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

# def search_for_fusion_information(request,c_line,algorithm,fusion_description,predicted_effect1,predicted_effect2):
#     response = {}
#     header = get_header()
#     
#     # recupero fusioni nella linea cellulare
#     fusions = []
#     
#     #linea cellulare non specificata, predicted effect specificato (0001)
#     if c_line == "ALL" and algorithm == "ALL" and fusion_description == "ALL" and predicted_effect1 != "ALL" and predicted_effect2 !="ALL":
#         for c_l in CellLine.nodes.all():
#             for fusion in c_l.happen:
#                 predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
#                 predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
#                 if predicted_effect_1 == predicted_effect1 and predicted_effect_2 == predicted_effect2:
#                     fusions.append(fusion)
#     #linea cellulare non specificata, mapping algorithm non specificato, fusion description specificata, predicted effect non specificato (0010)
#     elif c_line == "ALL" and algorithm == "ALL" and fusion_description != "ALL" and predicted_effect1 == "ALL" and predicted_effect2 =="ALL":
#         for c_l in CellLine.nodes.all():
#             for fusion in c_l.happen:
#                 if fusion_description in fusion.description:
#                     fusions.append(fusion)
#     #linea cellulare non specificata, mapping algorithm non specificato, fusion description specificata, predicted effect specificato (0011)
#     elif c_line == "ALL" and algorithm == "ALL" and fusion_description != "ALL" and predicted_effect1 != "ALL" and predicted_effect2 !="ALL":
#         for c_l in CellLine.nodes.all():
#             for fusion in c_l.happen:
#                 predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
#                 predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
#                 if fusion_description in fusion.description and predicted_effect_1 == predicted_effect1 and predicted_effect_2 == predicted_effect2:
#                     fusions.append(fusion)
#     #linea cellulare non specificata, mapping algorithm specificato, fusion description non specificata, predicted effect non speficicato (0100)
#     elif c_line == "ALL" and algorithm != "ALL" and fusion_description == "ALL" and predicted_effect1 == "ALL" and predicted_effect2 =="ALL":
#         for c_l in CellLine.nodes.all():
#             for fusion in c_l.happen:
#                 if algorithm in fusion.fusion_finding_method:
#                     fusion.append(fusion)
#     #linea cellulare non specificata, mapping algorithm specificato, fusion description non specificata, predicted effect specificato (0101)
#     elif c_line == "ALL" and algorithm != "ALL" and fusion_description == "ALL" and predicted_effect1 != "ALL" and predicted_effect2 !="ALL":
#         for c_l in CellLine.nodes.all():
#             for fusion in c_l.happen:
#                 predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
#                 predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
#                 if algorithm in fusion.fusion_finding_method and predicted_effect_1 == predicted_effect1 and predicted_effect_2 == predicted_effect2:
#                     fusion.append(fusion)
#     #linea cellulare non specificata, mapping algorithm specificato, fusion description specificata, predicted effect non specificato (0110)
#     elif c_line == "ALL" and algorithm != "ALL" and fusion_description != "ALL" and predicted_effect1 == "ALL" and predicted_effect2 =="ALL":
#         for c_l in CellLine.nodes.all():
#             for fusion in c_l.happen:
#                 if algorithm in fusion.fusion_finding_method and fusion_description in fusion.description:
#                     fusion.append(fusion)
#     #linea cellulare non specificata, mapping algorithm specificato, fusion description specificata, predicted effect specificato (0111)
#     elif c_line == "ALL" and algorithm != "ALL" and fusion_description != "ALL" and predicted_effect1 != "ALL" and predicted_effect2 !="ALL":
#         for c_l in CellLine.nodes.all():
#             for fusion in c_l.happen:
#                 predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
#                 predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
#                 if algorithm in fusion.fusion_finding_method and fusion_description in fusion.description and predicted_effect_1 == predicted_effect1 and predicted_effect_2 == predicted_effect2:
#                     fusion.append(fusion)
#     #linea cellulare specificata, il resto non specificato (1000) tutte le fusioni di una linea cellulare (analogo a search for cell_line)
#     elif c_line != "ALL" and algorithm == "ALL" and fusion_description == "ALL" and predicted_effect1 == "ALL" and predicted_effect2 =="ALL":
#         for fusion in CellLine.nodes.get(cell_line = c_line).happen:
#             fusions.append(fusion)     
#     #linea cellulare specificata, mapping algorithm non specificato, fusion description non specificata, predicted effect specificato (1001)
#     elif c_line != "ALL" and algorithm == "ALL" and fusion_description == "ALL" and predicted_effect1 != "ALL" and predicted_effect2 !="ALL":
#         for fusion in CellLine.nodes.get(cell_line = c_line).happen:
#             predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
#             predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
#             if predicted_effect_1 == predicted_effect1 and predicted_effect_2 == predicted_effect2:
#                 fusions.append(fusion)
#     #linea cellulare specificata, mapping algorithm non specificato, fusion description specificata, predicted effect non specificato (1010)
#     elif c_line != "ALL" and algorithm == "ALL" and fusion_description != "ALL" and predicted_effect1 == "ALL" and predicted_effect2 =="ALL":
#         for fusion in CellLine.nodes.get(cell_line = c_line).happen:
#             if fusion_description in fusion.description:
#                 fusions.append(fusion)
#     #linea cellulare specificata, mapping algorithm non specificato, fusion description specificata, predicted effect specificato (1011)
#     elif c_line != "ALL" and algorithm == "ALL" and fusion_description != "ALL" and predicted_effect1 != "ALL" and predicted_effect2 !="ALL":
#         for fusion in CellLine.nodes.get(cell_line = c_line).happen:
#                 predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
#                 predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
#                 if fusion_description in fusion.description and predicted_effect_1 == predicted_effect1 and predicted_effect_2 == predicted_effect2:
#                     fusions.append(fusion)
#     #linea cellulare specificata, mapping algorithm specificato, fusion description non specificato, predicted effect non specificato (1100)
#     elif c_line != "ALL" and algorithm != "ALL" and fusion_description == "ALL" and predicted_effect1 == "ALL" and predicted_effect2 =="ALL":
#         for fusion in CellLine.nodes.get(cell_line = c_line).happen:
#             if algorithm in fusion.fusion_finding_method:
#                     fusion.append(fusion)
#     #linea cellulare specificata, mapping algorithm specificato, fusion description non specificato, predicted effect specificato (1101)
#     elif c_line != "ALL" and algorithm != "ALL" and fusion_description == "ALL" and predicted_effect1 != "ALL" and predicted_effect2 !="ALL":
#         for fusion in CellLine.nodes.get(cell_line = c_line).happen:
#             predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
#             predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
#             if algorithm in fusion.fusion_finding_method and predicted_effect_1 == predicted_effect1 and predicted_effect_2 == predicted_effect2:
#                 fusion.append(fusion)
#     #linea cellulare specificata, mapping algorithm specificato, fusion description specificata, predicted effect non specificato (1110)
#     elif c_line != "ALL" and algorithm != "ALL" and fusion_description != "ALL" and predicted_effect1 == "ALL" and predicted_effect2 =="ALL":
#         for fusion in CellLine.nodes.get(cell_line = c_line).happen:
#             for fusion in c_l.happen:
#                 if algorithm in fusion.fusion_finding_method and fusion_description in fusion.description:
#                     fusion.append(fusion)
#     else:
#     #tutti i dati specificati (1111)
#         for fusion in CellLine.nodes.get(cell_line = c_line).happen:
#             predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
#             predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
#             if (algorithm in fusion.fusion_finding_method) and (fusion_description in fusion.description) and (predicted_effect1 == predicted_effect_1) and (predicted_effect2 == predicted_effect_2):
#                 fusions.append(fusion)
# 
#     rows = build_rows(fusions)
#     #print(rows)
# 
#     response['rows'] = {"header": header, "items": rows}
#     return HttpResponse(json.dumps(response))

def build_rows(fusions):
    rows = []
    
    for myfusion in fusions:
        # recupero dati cell line
        cellLine = myfusion.fromCellLineToFusion.all()[0].cell_line
        disease = ""
        acronym = ""
        ccle_infos = get_ccle_infos()
        for row in ccle_infos:
            if cellLine in row:
                disease = row[3]
                acronym = row[2]
        #recupero dati dai geni
        gene1 = myfusion.fromGeneToFusion.all()[0]
        gene2 = myfusion.with_gene.all()[0]
       
      
        #check
        fc_flag = "NO FC"
        es_flag = "NO ES"
        th_flag = "NO TH"
        
        if(myfusion.with_fc_script):
            fc_flag = '{ type: "button", action: "dialog", url: fusioncatcher/'+str(myfusion.fusion_id)+' }'
        if(myfusion.with_eric_script):
            es_flag = '{ type: "button", action: "dialog", url: ericscript/'+str(myfusion.fusion_id)+' }'
        if(myfusion.with_eric_script):
            th_flag = '{ type: "button", action: "dialog", url: tophat/'+str(myfusion.fusion_id)+' }'
        
        rows.append([disease,acronym,cellLine,gene1.symbol,gene2.symbol,fc_flag,es_flag])
        
        #print(gene1.symbol+" "+gene2.symbol)
        #print(fc_fusions)   
        #print(es_fusions)
        #print("\n\n\n")
        
    return rows


#BUILD_FC_TABLE E BUILD_ES_TABLE per costruire la tabella dei dettagli per ogni coppia. Quando premo il pulsante relativo all'algoritmo deve stamparmi tutte le fusioni relative

def build_fc_table(request,fus_id):
    
    response = {}
    header = get_fc_header()
    
    fusions = []
    
#     for fcFusion in Fusion.nodes.get(fusion_id = fus_id).with_fc_script:
#         fusions.append(fcFusion)
        
    for fusion in Fusion.nodes.all():
        if fusion.fusion_id == int(fus_id):
            if fusion.with_fc_script:
                fusions.append(fusion)
        
        
    rows = build_fc_rows(fusions)
    
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def build_es_table(request,fus_id):
    
    response = {}
    header = get_es_header()
    
    fusions = []
    
#     for esFusion in Fusion.nodes.get(fusion_id = fus_id).with_eric_script:
#         fusions.append(esFusion)

    for fusion in Fusion.nodes.all():
        if fusion.fusion_id == int(fus_id):
            if fusion.with_eric_script:
                fusions.append(fusion)
        
    rows = build_es_rows(fusions)
    
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def build_th_table(request,fus_id):
    
    response = {}
    header = get_tophat_header()
    
    fusions = []
    
#     for esFusion in Fusion.nodes.get(fusion_id = fus_id).with_eric_script:
#         fusions.append(esFusion)

    for fusion in Fusion.nodes.all():
        if fusion.fusion_id == int(fus_id):
            if fusion.with_tophat_script:
                fusions.append(fusion)
        
    rows = build_tophat_rows(fusions)
    
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def build_fc_rows(fusions):
    
    rows = []
    # ora che ho solo le fusioni FC interessate recupero le informazioni e mi costruisco la riga
    fc_fusions = []
    
    for fus in fusions:
        for fcfusion in fus.with_fc_script:
            fc_fusions.append(fcfusion)
    for myfusion in fc_fusions:
        gene1 = myfusion.fromFusionToFusionCatcher.all()[0].fromGeneToFusion.all()[0]
        gene2 = myfusion.fromFusionToFusionCatcher.all()[0].with_gene.all()[0]
        
        strand_1 = myfusion.strand_1
        predicted_effect_1 = myfusion.predicted_effect_1
        strand_2 = myfusion.strand_2
        predicted_effect_2 = myfusion.predicted_effect_2
        
        #recupero cromosomi 
        chromosome1 = gene1.fromChromosomeToGene.all()[0]
        chromosome2 = gene2.fromChromosomeToGene.all()[0]
        fusion_point_1 = myfusion.fusion_point_1
        fusion_point_2 = myfusion.fusion_point_2

        #recupero esoni        
        exon1 = []
        exon2 = []

        for exon in myfusion.at_exon:
            if exon.in_gene.filter(symbol__exact=gene1.symbol):
                exon1 = exon
            if exon.in_gene.filter(symbol__exact=gene2.symbol):
                exon2 = exon
        #recupero trascritti e proteine
        transcript_couples = []
        proteins = []
        for couple in myfusion.with_trans_couple:
            transcript1 = couple.fromTranscriptToCouple.all()[0]
            transcript1_position = couple.fromTranscriptToCouple.relationship(transcript1).position
            transcript2 = couple.with_other_transcript.all()[0]
            transcript2_position = couple.with_other_transcript.relationship(transcript2).position
            
            transcript_couples.append(transcript1.transcript+":"+str(transcript1_position)+" - "+transcript2.transcript+":"+str(transcript2_position))
            if couple.with_protein.all():
                proteins.append(couple.with_protein.all()[0].protein)
            
        #costruisco la riga
        row = []
        row.append(gene1.symbol+" - "+gene2.symbol)
        if exon1 or exon2:
            row.append(exon1.exon+" - "+exon2.exon)
        else:
            row.append("No exons")
        row.append([chromosome1.chromosome+":"+str(fusion_point_1)+":"+strand_1, chromosome2.chromosome+":"+str(fusion_point_2)+":"+strand_2])
        row.append(myfusion.description)
        row.append(myfusion.common_mapping_reads)
        row.append(myfusion.spanning_pairs)
        row.append(myfusion.spanning_unique_reads)
        row.append(myfusion.longest_anchor_found)
        row.append(myfusion.fusion_finding_method)
        row.append(myfusion.fusion_sequence)
        if(predicted_effect_1!=predicted_effect_2):
            row.append(predicted_effect_1+"/"+predicted_effect_2)
        else:
            row.append(predicted_effect_1)
        row.append(transcript_couples)
        row.append(proteins)   
        #print(row)
        rows.append(row)
    return rows

#SCRIVERE BUILD_ES_ROWS
def build_es_rows(fusions):
    rows = []
    #recupero dati degli eventi di fusione
        
    es_fusions = []
    
    for fus in fusions:
        for esfusion in fus.with_eric_script:
            es_fusions.append(esfusion)
    
    for myfusion in es_fusions:
        gene1 = myfusion.fromFusionToEricScript.all()[0].fromGeneToFusion.all()[0]
        gene2 = myfusion.fromFusionToEricScript.all()[0].with_gene.all()[0]
        
        #recupero cromosomi 
        chromosome1 = gene1.fromChromosomeToGene.all()[0]
        chromosome2 = gene2.fromChromosomeToGene.all()[0]    
        breakpoint1 = myfusion.breakpoint_1
        breakpoint2 = myfusion.breakpoint_2
        strand1 = myfusion.strand_1
        strand2 = myfusion.strand_2
        
        crossing_reads = myfusion.crossing_reads
        spanning_reads = myfusion.spanning_reads
        mean_intersize = myfusion.mean_intersize
        homology = myfusion.homology
        fusion_type = myfusion.fusion_type
        junction_sequence = myfusion.junction_sequence
        gene_expr_1 = myfusion.gene_expr_1
        gene_expr_2 = myfusion.gene_expr_2
        gene_expr_fused = myfusion.gene_expr_fused
        es = myfusion.es
        gjs = myfusion.gjs
        us = myfusion.us
        eric_score = myfusion.eric_score
        
        #costruisco la riga
        row = []
        row.append(gene1.symbol+" - "+gene2.symbol)
        row.append([chromosome1.chromosome+":"+str(breakpoint1)+":"+strand1, chromosome2.chromosome+":"+str(breakpoint2)+":"+strand2])
        row.append(crossing_reads)
        row.append(spanning_reads)
        row.append(mean_intersize)
        row.append(homology)
        row.append(fusion_type)
        row.append(junction_sequence)
        row.append(gene_expr_1)
        row.append(gene_expr_2)
        row.append(gene_expr_fused)
        row.append(es)
        row.append(gjs)
        row.append(us)
        row.append(eric_score)
        
        rows.append(row)

    return rows

def build_tophat_rows(fusions):
    rows = []
    #recupero dati degli eventi di fusione
        
    th_fusions = []
    
    for fus in fusions:
        for thfusion in fus.with_tophatscript:
            th_fusions.append(thfusion)
    
    for myfusion in th_fusions:
        gene1 = myfusion.fromFusionToTophatScript.all()[0].fromGeneToFusion.all()[0]
        gene2 = myfusion.fromFusionToTophatScript.all()[0].with_gene.all()[0]
        
        #recupero cromosomi 
        chromosome1 = gene1.fromChromosomeToGene.all()[0]
        chromosome2 = gene2.fromChromosomeToGene.all()[0]    
        left_coord = myfusion.left_coord
        right_coord = myfusion.right_coord
        
        spanning_reads = myfusion.spanning_reads
        spanning_mate_pairs = myfusion.spanning_mate_pairs
        spanning_mate_pairs_end = myfusion.spanning_mate_pairs_end
        
        #costruisco la riga
        row = []
        row.append(gene1.symbol+" - "+gene2.symbol)
        row.append([chromosome1.chromosome+":"+str(left_coord), chromosome2.chromosome+":"+str(right_coord)])
        row.append(spanning_reads)
        row.append(spanning_mate_pairs)
        row.append(spanning_mate_pairs_end)

        
        rows.append(row)

    return rows
    

def search_viruses(request,c_line,vir):
    response = {}
    rows = []
    header = ["Cell line",
        "Virus/bacteria/phage name",
        "GI",
        "NC"]
    
    #recupero virus nella linea cellulare
    #se ho specificato solo il virus, cerca per tutte le linee cellulari
    if c_line == "ALL" and vir != "":
        for c_l in CellLine.nodes.all():
            for virus in c_l.with_viruses:
                if vir in virus.name or vir == virus.gi or vir == virus.ref:
                    rows.append([c_l.cell_line, virus.name, virus.gi, virus.ref])
    #se ho specificato solo la linea cellulare, cerca tutti i virus per la linea cellulare
    elif c_line !="" and vir == "trolo":
        for virus in CellLine.nodes.get(cell_line = c_line).with_viruses:
            rows.append([c_line, virus.name, virus.gi, virus.ref])
    #se ho specificato sia virus che linea cellulare, cerca quel virus per quella linea cellulare
    elif c_line != "" and vir != "":
        for virus in CellLine.nodes.get(cell_line = c_line).with_viruses:
            if vir in virus.name or vir == virus.gi or vir == virus.ref:
                rows.append([c_line, virus.name, virus.gi, virus.ref])
                
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))
        
def generate_statistics(request):
    cellLine_gene = {(('CellLine','cell_line'),('Gene','symbol')):[[CellLine,'HAPPEN',Fusion,'HAD',Gene],[CellLine,'HAPPEN',Fusion,'WITH',Gene]]}
    cellLine_chromosome = {(('CellLine','cell_line'),('Chromosome','chromosome')):[[CellLine,'HAPPEN',Fusion,'AT_CHROMOSOME',Chromosome]]}
    cellLine_fusion = {(('CellLine','cell_line'),('Fusion','fusion_id')):[[CellLine,'HAPPEN',Fusion]]}
    chromosome_cellLine = {(('Chromosome','chromosome'),('CellLine','cell_line')):[[Chromosome,'AT_CHROMOSOME',Fusion,'HAPPEN',CellLine]]}
    chromosome_fusion = {(('Chromosome','chromosome'),('Fusion','fusion_id')):[[Chromosome,'AT_CHROMOSOME',Fusion]]}
    chromosome_gene = {(('Chromosome','chromosome'),('Gene','symbol')):[[Chromosome,'OF_GENE',Gene]]}
    gene_cellLine = {(('Gene','symbol'),('CellLine','cell_line')):[[Gene,'HAD',Fusion,'HAPPEN',CellLine],[Gene,'WITH',Fusion,'HAPPEN',CellLine]]}
    gene_chromosome = {(('Gene','symbol'),('Chromosome','chromosome')):[[Gene,'AT_CHROMOSOME',Chromosome]]}
    gene_fusion = {(('Gene','symbol'),('Fusion','fusion_id')):[[Gene,'HAD',Fusion],[Gene,'WITH',Fusion]]}
    
    stats = [cellLine_gene, cellLine_chromosome, cellLine_fusion, chromosome_cellLine,chromosome_fusion, chromosome_gene, gene_cellLine, gene_chromosome, gene_fusion]
    for pair in stats:
        gen_stat_file(pair)
    #gen_stat_file(gene_cellLine)
    
    return HttpResponse()    
    
def gen_stat_file(pairs):
    ids = list(pairs.keys())[0]
    paths = list(pairs.values())
    paths = paths[0]
    #print(ids)
    #print(paths)
    startnode = paths[0][0]
    #print(startnode)
    file =  open(ids[0][0]+'_'+ids[1][0]+'.csv','w')
    writer = csv.writer(file, lineterminator='\n')
    writer.writerow([ids[0][0],ids[1][0]])
    for x in startnode.nodes.all():
        print(x)
        nodes = {}
        for path in paths:
            #print(path)
            #print(len(path))
            query = "match (" + path[0].__name__ + "{"+ids[0][1]+":'"+eval("x."+ids[0][1])+"'})-"
            for i in range(1,len(path)-2,2):
                query = query + "[:" +  path[i] + "]-(" + path[i+1].__name__ + ")-" 
                #print(path[i+1])
            query = query + "[:" + path[len(path)-2] + "]-(a:" + path[(len(path)-1)].__name__+") return distinct a"
            #print(query)
            if(db.cypher_query(query)[0]):
                for row in db.cypher_query(query)[0]:
                    #print(row[0].properties[eval("'"+ids[1][1]+"'")])
                    if row[0].properties[eval("'"+ids[1][1]+"'")] not in nodes:
                        #print(row[1].properties[eval("'"+ids[1][1]+"'")])
                        nodes[row[0].properties[eval("'"+ids[1][1]+"'")]] = row[0].properties[eval("'"+ids[1][1]+"'")]
        #print(len(nodes))
        writer.writerow([eval("x."+ids[0][1]),len(nodes)])
    file.close()
    
def get_ccle_infos():
    header = ["ID","Cell Line","Disease","Disease name"]
    rows = []
    response = {}
    
    txt_file = open("ccle_ids.txt", "r")
    next(txt_file)
    for line in txt_file:
        words = line.split("\t")
        rows.append([words[0].replace(" ",""),words[1],words[2],words[3].replace("\n","")])
        #print(words)
    
    #response['rows'] = {"header": header, "items": rows}

    return rows

def get_header():
    return ["Cancer",
            "Acronym",
            "CCLE",
            "Gene 1",
            "Gene 2",
            "Fusion Catcher",
            "Ericscript",
            "Tophat"]
    
def get_fc_header():
    return ["Gene pair symbols",
        "Exon pair",
        "Chromosome : fusion point : strand",
        "Description",
        "Counts of common mapping reads",
        "Spanning pairs",
        "Spanning unique reads",
        "Longest anchor found",
        "Fusion finding method",
        "Fusion sequence",
        "Predicted effect",
        "Predicted fused transcripts",
        "Predicted fused proteins"]
    
def get_es_header():
    return ["Gene pair symbols",
        "Chromosome : breakpoint : strand",
        "Crossing reads",
        "Spanning reads",
        "Mean intersize",
        "Homology",
        "Fusion type",
        "Junction sequence",
        "Gene expr 1",
        "Gene expr 2",
        "Gene expr fused",
        "Es",
        "Gjs",
        "Us",
        "Eric score"]
    
def get_tophat_header():
    return ["Gene pair symbols",
            "Chromosome : coordinate",
            "Spanning reads",
            "Spanning mate pairs",
            "Spanning mate pairs where one end spans a fusion"]
    
    #prendo in input una stringa che e' il nome della malattia, mi ricavo le linee cellulari corrispondenti e mi ricavo la tabella relativa
def get_cell_line_from_disease(disease):
    ccle_infos = get_ccle_infos()
    cls = []
    for row in ccle_infos:
        if disease in row:
            cls.append(row[0])
    return cls

def search_for_disease(request,disease):
    cls = get_cell_line_from_disease(disease)
    fusions = []
    for cl in cls:
        response = {}
        header = get_header()
        for fusion in CellLine.nodes.get(cell_line = cl).happen:
            fusions.append(fusion)
        
    rows = build_rows(fusions)
    
    response['rows'] = {"header": header, "items": rows}    
        
    return HttpResponse(json.dumps(response))    
        
    