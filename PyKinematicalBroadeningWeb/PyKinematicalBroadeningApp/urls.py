from django.urls import path
from . import views

urlpatterns = [
    path('admin/', admin.site.urls),
    path('broaden_spectrum/', views.broaden_spectrum, name='broaden_spectrum'),
]
