"""
Backend girişi: config.yaml'ı okuyup iSAM2 artımlı BA çalıştırır; run_frontend.py dump'ını
tüketir, trajectory'yi kestirir, GT'ye ankorlar, CSV + ATE yazar.

    python run_backend.py                        # config.yaml'daki ayarlarla
    python run_backend.py --stereo_sigma 0.35    # tek seferlik CLI override
"""
from gtsam_backend import main

if __name__ == "__main__":
    main()
