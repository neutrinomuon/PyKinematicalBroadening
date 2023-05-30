from django.contrib import admin
from django.urls import path
from PyKinematicalBroadeningApp.views import broaden_spectrum

urlpatterns = [
    path('admin/', admin.site.urls),
    path('broaden_spectrum/', broaden_spectrum, name='broaden_spectrum')
]


