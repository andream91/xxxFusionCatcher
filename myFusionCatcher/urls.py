"""myFusionCatcher URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.10/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import url
from django.contrib import admin
from app import views

urlpatterns = [
    url(r'^admin/', admin.site.urls),
    url(r'^chromosomes/(.+)/(.+)/(.+)/', views.search_for_chromosome, name='search_for_chromosome'),
    url(r'^cell_line/(.+)/', views.search_for_cell_line, name='search_for_cell_line'),
    url(r'^gene/(.+)/(.+)/(.+)/', views.search_for_gene, name='search_for_gene'),
    url(r'^exon/(.+)/(.+)/(.+)/', views.search_for_exon, name='search_for_exon'),
    url(r'^transcript/(.+)/(.+)/(.+)/', views.search_for_transcript, name='search_for_transcript'),
    #url(r'^fusion_information/(.+)/(.+)/(.+)/(.+)/(.+)/', views.search_for_fusion_information, name='search_for_fusion_information'),
    url(r'^disease/(.+)/', views.search_for_disease, name='search_for_disease'),
    url(r'^virus/(.+)/(.+)/', views.search_viruses, name='search_viruses'),
    url(r'^genstats/', views.generate_statistics, name='generate_statistics'),
    url(r'^fusioncatcher/(.+)/', views.build_fc_table, name='build_fc_table'),
    url(r'^ericscript/(.+)/', views.build_es_table, name='build_es_table'),
    url(r'^tophat/(.+)/', views.build_th_table, name='build_th_table'),

]