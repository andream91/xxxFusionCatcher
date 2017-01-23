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
    url(r'^chromosomes/(.+)/(.+)/(.+)/(.+)/', views.search_for_chromosome, name='search_for_chromosome'),
    url(r'^cell_line/(.+)/', views.search_for_cell_line, name='search_for_cell_line'),
    url(r'^gene/single/(.+)/(.+)/', views.search_for_single_gene, name='search_for_single_gene'),
    url(r'^gene/pair/(.+)/(.+)/(.+)/', views.search_for_pair_gene, name='search_for_pair_gene'),
    url(r'^exon/single/(.+)/(.+)/', views.search_for_single_exon, name='search_for_single_exon'),
    url(r'^exon/pair/(.+)/(.+)/(.+)/', views.search_for_pair_exon, name='search_for_pair_exon'),
    url(r'^transcript/single/(.+)/(.+)/', views.search_for_single_transcript, name='search_for_single_transcript'),
    url(r'^transcript/pair/(.+)/(.+)/(.+)/', views.search_for_pair_transcript, name='search_for_pair_transcript'),
    url(r'^fusion_information/(.+)/(.+)/(.+)/(.+)/(.+)/', views.search_for_fusion_information, name='search_for_fusion_information'),
]
