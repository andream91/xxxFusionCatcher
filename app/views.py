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

def search_for_chromosome(request,c_line,chromos):
    response = {}
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    #se ho specificato solo il cromosoma, cerco tutte le fusioni in tutte le linee cellulari che coinvolgono il cromosoma
    if c_line == "ALL" and chromos != "":
        c = Chromosome.nodes.get(chromosome = chromos)
        #for fcfusion in c.fromFCFusiontoChromosome:
        #    for fusion in fcfusion.fromFusionToFusionCatcher:
        #        fusions.append(fusion)
        #print(c.fromFCFusiontoChromosome[0])

        print(c.fromESFusiontoChromosome)
        
        #for esfusion in c.fromESFusiontoChromosome:
        #    print("BOH")
            #for fusion in esfusion.fromFusionToEricScript:
            #    fusions.append(fusion)
    #se ho specificato solo la linea cellulare, cerco tutte le fusioni che avvengono in uqella linea cellulare nel determinato intervallo
    #elif c_line != "" and chromos == "":
    #    for fusion in CellLine.nodes.get(cell_line = c_line).happen:
    #            fusions.append(fusion)
    #se ho sia linea cellulare che cromosoma specificati, cerco tutte le fusioni nella linea cellulare che coinvolgono il cromosoma
    #elif c_line != "" and chromos != "":
    #    for fusion in CellLine.nodes.get(cell_line = c_line).happen:
    #            if fusion.at_chromosome.filter(chromosome__exact=chromos):
    #                fusions.append(fusion)
    
    rows = build_rows(fusions)
    #print(rows)
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))


def OLD_search_for_chromosome(request,c_line,chromos,start_point,end_point):
    response = {}
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    #se ho specificato solo il cromosoma, cerco tutte le fusioni in tutte le linee cellulari che coinvolgono il cromosoma
    if c_line == "ALL" and chromos != "":
        c = Chromosome.nodes.get(chromosome = chromos)
        for fusion in c.fromFusiontoChromosome:
            if fusion.fusion_point_1 >= start_point and fusion.fusion_point_2 <= end_point:
                fusions.append(fusion)
    #se ho specificato solo la linea cellulare, cerco tutte le fusioni che avvengono in uqella linea cellulare nel determinato intervallo
    elif c_line != "" and chromos == "":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion.fusion_point_1 >= start_point and fusion.fusion_point_2 <= end_point:
                fusions.append(fusion)
    #se ho sia linea cellulare che cromosoma specificati, cerco tutte le fusioni nella linea cellulare che coinvolgono il cromosoma
    elif c_line != "" and chromos != "":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion.fusion_point_1 >= int(start_point) and fusion.fusion_point_2 <= int(end_point):
                if fusion.at_chromosome.filter(chromosome__exact=chromos):
                    fusions.append(fusion)
    
    rows = build_rows(fusions)
    #print(rows)
    
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_gene(request,c_line,gene_one,gene_two):
    response = {}
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    
    #se ho specificato il gene ma non la linea cellulare, mostro tutte le fusioni che coinvolgono il gene
    if c_line == "ALL" and gene_one != "" and gene_two == "":
        if "ENSG" in gene_one:
            g = Gene.nodes.get(gene_id = gene_one)
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
    elif c_line == "ALL" and gene_one != "" and gene_two != "":
        if "ENSG" in gene_one and "ENSG" in gene_two:
            g = Gene.nodes.get(gene_id = gene_one)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].gene_id == gene_two:
                        fusions.append(fusion)
        
            g = Gene.nodes.get(gene_id = gene_two)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].gene_id == gene_one:
                        fusions.append(fusion)
        elif "ENSG" in gene_one and "ENSG" not in gene_two:
            g = Gene.nodes.get(gene_id = gene_one)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].symbol == gene_two:
                        fusions.append(fusion)
        
            g = Gene.nodes.get(symbol = gene_two)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].gene_id == gene_one:
                        fusions.append(fusion)
        elif "ENSG" not in gene_one and "ENSG" in gene_two:
            g = Gene.nodes.get(symbol = gene_one)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].gene_id == gene_two:
                        fusions.append(fusion)
        
            g = Gene.nodes.get(gene_id = gene_two)
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
    elif c_line != "ALL" and gene_one == "" and gene_two == "":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            fusions.append(fusion)     
    #se ho specificato sia la linea cellulare che il gene, cerco tutte le fusioni che coinvolgono quel gene nella linea cellulare
    elif c_line != "ALL" and gene_one != "" and gene_two == "":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion.fromGeneToFusion.filter(symbol__exact=gene_one) or fusion.fromGeneToFusion.filter(gene_id__exact=gene_one) or fusion.with_gene.filter(symbol__exact=gene_one) or fusion.with_gene.filter(gene_id__exact=gene_one):
                fusions.append(fusion)
    elif c_line != "ALL" and gene_one != "" and gene_two != "":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if (fusion.fromGeneToFusion.filter(symbol__exact=gene_one) and fusion.with_gene.filter(symbol__exact=gene_two)) or (fusion.fromGeneToFusion.filter(symbol__exact=gene_two) and fusion.with_gene.filter(symbol__exact=gene_one)) or (fusion.fromGeneToFusion.filter(gene_id__exact=gene_one) and fusion.with_gene.filter(gene_id__exact=gene_two)) or (fusion.fromGeneToFusion.filter(gene_id__exact=gene_two) and fusion.with_gene.filter(gene_id__exact=gene_one)):
                fusions.append(fusion)
            
    rows = build_rows(fusions)
    #print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))



def search_for_exon(request,c_line,exon_one,exon_two):
    response = {}
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    
    #cell_line all, primo esone si, secondo esone no, tutte le fusioni che coinvolgono questo esone
    if c_line == "ALL" and exon_one!="" and exon_two == "":
        e = Exon.nodes.get(exon = exon_one)
        for fusion in e.fromFusionToExon:
            fusions.append(fusion)
    #cell_line all, primo esone si, secondo esone si (011), tutte le fusioni che convolgono la coppia di esoni
    elif c_line == "ALL" and exon_one!="" and exon_two!="":
        e = Exon.nodes.get(exon = exon_one)
        for fusion in e.fromFusionToExon:
            if fusion.at_exon.filter(exon__exact=exon_two):
                fusions.append(fusion)
    #cell_line si, primo esone no, secondo esone no, ANALOGO A SEARCH FOR CELL_LINE
    elif c_line!="ALL" and exon_one=="" and exon_two=="":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            fusions.append(fusion)     
    #cell_line si, primo esone si, secondo esone no 
    elif c_line != "ALL" and exon_one != "" and exon_two == "":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion.at_exon.filter(exon__exact=exon_one):
                fusions.append(fusion)
    #cell_line si, primo esone si, secondo esone si 
    elif c_line != "ALL" and exon_one!="" and exon_two != "":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion.at_exon.filter(exon__exact=exon_one) and fusion.at_exon.filter(exon__exact=exon_two):
                    fusions.append(fusion)

    rows = build_rows(fusions)
    #print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_transcript(request,c_line,transcript_one,transcript_two):
    response = {}
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    
    #linea cellulare non specificata, un trascritto specificato
    if c_line=="ALL" and transcript_one!="" and transcript_two=="":
        for couple in Transcript.nodes.get(transcript = transcript_one).fromCoupleToTranscript:
            for fusion in couple.fromFusionToCouple:
                fusions.append(fusion)
    #linea cellulare non specificata, coppia di trascritti specificata
    elif c_line=="ALL" and transcript_one!="" and transcript_two!="":
        for couple in Transcript.nodes.get(transcript = transcript_one).fromCoupleToTranscript:
            if couple.with_other_transcript.filter(transcript__exact=transcript_two) or couple.fromTranscriptToCouple.filter(transcript__exact=transcript_two):
                for fusion in couple.fromFusionToCouple:
                    fusions.append(fusion)
        for couple in Transcript.nodes.get(transcript = transcript_two).fromCoupleToTranscript:
            if couple.with_other_transcript.filter(transcript__exact=transcript_one) or couple.fromTranscriptToCouple.filter(transcript__exact=transcript_one):
                for fusion in couple.fromFusionToCouple:
                    fusions.append(fusion)
    #linea cellulare specificata, un trascritto specificato:
    elif c_line!="ALL" and transcript_one!="" and transcript_two=="trolo":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            for couple in fusion.with_trans_couple:
                if couple.fromTranscriptToCouple.filter(transcript__exact=transcript_one) or couple.with_other_transcript.filter(transcript__exact=transcript_one):
                    fusions.append(fusion)
    #linea cellulare specificata, coppia di trascritti specificata
    elif c_line!="ALL" and transcript_one!="" and transcript_two!="":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            for couple in fusion.with_trans_couple:
                if (couple.fromTranscriptToCouple.filter(transcript__exact=transcript_one) and couple.with_other_transcript.filter(transcript__exact=transcript_two)) or (couple.fromTranscriptToCouple.filter(transcript__exact=transcript_two) and couple.with_other_transcript.filter(transcript__exact=transcript_one)):
                    fusions.append(fusion)
        

    rows = build_rows(fusions)
    #print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_fusion_information(request,c_line,algorithm,fusion_description,predicted_effect1,predicted_effect2):
    response = {}
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    
    #linea cellulare non specificata, predicted effect specificato (0001)
    if c_line == "ALL" and algorithm == "ALL" and fusion_description == "ALL" and predicted_effect1 != "ALL" and predicted_effect2 !="ALL":
        for c_l in CellLine.nodes.all():
            for fusion in c_l.happen:
                predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
                predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
                if predicted_effect_1 == predicted_effect1 and predicted_effect_2 == predicted_effect2:
                    fusions.append(fusion)
    #linea cellulare non specificata, mapping algorithm non specificato, fusion description specificata, predicted effect non specificato (0010)
    elif c_line == "ALL" and algorithm == "ALL" and fusion_description != "ALL" and predicted_effect1 == "ALL" and predicted_effect2 =="ALL":
        for c_l in CellLine.nodes.all():
            for fusion in c_l.happen:
                if fusion_description in fusion.description:
                    fusions.append(fusion)
    #linea cellulare non specificata, mapping algorithm non specificato, fusion description specificata, predicted effect specificato (0011)
    elif c_line == "ALL" and algorithm == "ALL" and fusion_description != "ALL" and predicted_effect1 != "ALL" and predicted_effect2 !="ALL":
        for c_l in CellLine.nodes.all():
            for fusion in c_l.happen:
                predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
                predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
                if fusion_description in fusion.description and predicted_effect_1 == predicted_effect1 and predicted_effect_2 == predicted_effect2:
                    fusions.append(fusion)
    #linea cellulare non specificata, mapping algorithm specificato, fusion description non specificata, predicted effect non speficicato (0100)
    elif c_line == "ALL" and algorithm != "ALL" and fusion_description == "ALL" and predicted_effect1 == "ALL" and predicted_effect2 =="ALL":
        for c_l in CellLine.nodes.all():
            for fusion in c_l.happen:
                if algorithm in fusion.fusion_finding_method:
                    fusion.append(fusion)
    #linea cellulare non specificata, mapping algorithm specificato, fusion description non specificata, predicted effect specificato (0101)
    elif c_line == "ALL" and algorithm != "ALL" and fusion_description == "ALL" and predicted_effect1 != "ALL" and predicted_effect2 !="ALL":
        for c_l in CellLine.nodes.all():
            for fusion in c_l.happen:
                predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
                predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
                if algorithm in fusion.fusion_finding_method and predicted_effect_1 == predicted_effect1 and predicted_effect_2 == predicted_effect2:
                    fusion.append(fusion)
    #linea cellulare non specificata, mapping algorithm specificato, fusion description specificata, predicted effect non specificato (0110)
    elif c_line == "ALL" and algorithm != "ALL" and fusion_description != "ALL" and predicted_effect1 == "ALL" and predicted_effect2 =="ALL":
        for c_l in CellLine.nodes.all():
            for fusion in c_l.happen:
                if algorithm in fusion.fusion_finding_method and fusion_description in fusion.description:
                    fusion.append(fusion)
    #linea cellulare non specificata, mapping algorithm specificato, fusion description specificata, predicted effect specificato (0111)
    elif c_line == "ALL" and algorithm != "ALL" and fusion_description != "ALL" and predicted_effect1 != "ALL" and predicted_effect2 !="ALL":
        for c_l in CellLine.nodes.all():
            for fusion in c_l.happen:
                predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
                predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
                if algorithm in fusion.fusion_finding_method and fusion_description in fusion.description and predicted_effect_1 == predicted_effect1 and predicted_effect_2 == predicted_effect2:
                    fusion.append(fusion)
    #linea cellulare specificata, il resto non specificato (1000) tutte le fusioni di una linea cellulare (analogo a search for cell_line)
    elif c_line != "ALL" and algorithm == "ALL" and fusion_description == "ALL" and predicted_effect1 == "ALL" and predicted_effect2 =="ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            fusions.append(fusion)     
    #linea cellulare specificata, mapping algorithm non specificato, fusion description non specificata, predicted effect specificato (1001)
    elif c_line != "ALL" and algorithm == "ALL" and fusion_description == "ALL" and predicted_effect1 != "ALL" and predicted_effect2 !="ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
            predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
            if predicted_effect_1 == predicted_effect1 and predicted_effect_2 == predicted_effect2:
                fusions.append(fusion)
    #linea cellulare specificata, mapping algorithm non specificato, fusion description specificata, predicted effect non specificato (1010)
    elif c_line != "ALL" and algorithm == "ALL" and fusion_description != "ALL" and predicted_effect1 == "ALL" and predicted_effect2 =="ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion_description in fusion.description:
                fusions.append(fusion)
    #linea cellulare specificata, mapping algorithm non specificato, fusion description specificata, predicted effect specificato (1011)
    elif c_line != "ALL" and algorithm == "ALL" and fusion_description != "ALL" and predicted_effect1 != "ALL" and predicted_effect2 !="ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
                predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
                predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
                if fusion_description in fusion.description and predicted_effect_1 == predicted_effect1 and predicted_effect_2 == predicted_effect2:
                    fusions.append(fusion)
    #linea cellulare specificata, mapping algorithm specificato, fusion description non specificato, predicted effect non specificato (1100)
    elif c_line != "ALL" and algorithm != "ALL" and fusion_description == "ALL" and predicted_effect1 == "ALL" and predicted_effect2 =="ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if algorithm in fusion.fusion_finding_method:
                    fusion.append(fusion)
    #linea cellulare specificata, mapping algorithm specificato, fusion description non specificato, predicted effect specificato (1101)
    elif c_line != "ALL" and algorithm != "ALL" and fusion_description == "ALL" and predicted_effect1 != "ALL" and predicted_effect2 !="ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
            predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
            if algorithm in fusion.fusion_finding_method and predicted_effect_1 == predicted_effect1 and predicted_effect_2 == predicted_effect2:
                fusion.append(fusion)
    #linea cellulare specificata, mapping algorithm specificato, fusion description specificata, predicted effect non specificato (1110)
    elif c_line != "ALL" and algorithm != "ALL" and fusion_description != "ALL" and predicted_effect1 == "ALL" and predicted_effect2 =="ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            for fusion in c_l.happen:
                if algorithm in fusion.fusion_finding_method and fusion_description in fusion.description:
                    fusion.append(fusion)
    else:
    #tutti i dati specificati (1111)
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
            predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
            if (algorithm in fusion.fusion_finding_method) and (fusion_description in fusion.description) and (predicted_effect1 == predicted_effect_1) and (predicted_effect2 == predicted_effect_2):
                fusions.append(fusion)

    rows = build_rows(fusions)
    #print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

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
       
        rows.append([disease,acronym,cellLine,gene1.symbol,gene2.symbol])
       
    #print(rows)   
    return rows



def build_fc_rows(fusions):
    
    rows = []
    # ora che ho solo le fusioni interessate recupero le informazioni e mi costruisco la riga
    for myfusion in fusions:
        # recupero cell line
        cellLine = myfusion.fromCellLineToFusion.all()[0].cell_line
        
        #recupero dati dai geni
        gene1 = myfusion.fromGeneToFusion.all()[0]
        strand_1 = myfusion.fromGeneToFusion.relationship(gene1).strand
        predicted_effect_1 = myfusion.fromGeneToFusion.relationship(gene1).predicted_effect
        gene2 = myfusion.with_gene.all()[0]
        strand_2 = myfusion.with_gene.relationship(gene2).strand
        predicted_effect_2 = myfusion.with_gene.relationship(gene2).predicted_effect

        #recupero cromosomi 
        chromosome1 = myfusion.at_chromosome.match(fusion_partner__exact="5'end").all()[0]
        chromosome2 = myfusion.at_chromosome.match(fusion_partner__exact="3'end").all()[0]
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
        row.append(cellLine)
        row.append(gene1.symbol+" - "+gene2.symbol)
        row.append(gene1.gene_id+" - "+gene2.gene_id)
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
        
def old_generate_statistics(request):
    pairs = {CellLine:('CellLine','cell_line'),Fusion:('Fusion','fusion_id'),Chromosome:('Chromosome','chromosome'),Gene:('Gene','symbol') }
    for node1,node_data1 in pairs.items():
        for node2,node_data2 in pairs.items():
            if node1 != node2:
                if node_data1[0]=="Fusion":
                    print("Fusion -no.") #fusion-cell_line è inutile e ci mette tempo, abolirei anche tutti quelli che iniziano con fusion
                else:
                    print(node_data1[0],node_data2[0])
                    file =  open(node_data1[0]+'_'+node_data2[0]+'.csv','w')
                    writer = csv.writer(file, lineterminator='\n')
                    writer.writerow([node_data1[0],node_data2[0]])
                    for x in node1.nodes.all():
                        #query = "match (x:"+node_data1[0])-[*..2]-(y:"+str(node_data2[0])+") return x, count(distinct y)"
                        query = "match (x:"+node_data1[0]+"{"+node_data1[1]+":'"+str(eval("x."+str(node_data1[1])))+"'})-[*..1]-(y:"+str(node_data2[0])+") return x, count(distinct y)"
                        if(db.cypher_query(query)[0]): #ho la linea cellulare vuota, machecazz?
                            #print([db.cypher_query(query)[0][0][0].properties[eval("'"+node_data1[1]+"'")],db.cypher_query(query)[0][0][1]])
                            writer.writerow([db.cypher_query(query)[0][0][0].properties[eval("'"+node_data1[1]+"'")],db.cypher_query(query)[0][0][1]])
                    file.close() 
        
    return HttpResponse()
    
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
            "Ericscript"]
    
def get_fc_header():
    return ["Cell line",
        "Gene pair symbols",
        "Gene pair EnsIDs",
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
    
    #prendo in input una stringa che è il nome della malattia, mi ricavo le linee cellulari corrispondenti e mi ricavo la tabella relativa
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
        
    