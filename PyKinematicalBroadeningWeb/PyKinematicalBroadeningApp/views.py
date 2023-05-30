from django.shortcuts import render
from .utils import broadening_gh

def broadening_view(request):
    if request.method == 'POST':
        # Get form data
        lambda_o = [float(x) for x in request.POST.get('lambda_o').split(',')]
        fluxes_o = [float(x) for x in request.POST.get('fluxes_o').split(',')]
        lambda_s = [float(x) for x in request.POST.get('lambda_s').split(',')]
        vc0_gals = float(request.POST.get('vc0_gals'))
        vd_sigma = float(request.POST.get('vd_sigma'))

        # Call broadening function
        fluxes_s, IsKeepOn = broadening_gh(lambda_o, fluxes_o, lambda_s, vc0_gals, vd_sigma)

        # Render results
        return render(request, 'results.html', {'fluxes_s': fluxes_s, 'IsKeepOn': IsKeepOn})
    else:
        # Render form
        return render(request, 'form.html')
