"""
turb — Türbülans-Farkında Görsel-Ataletsel Odometri araç takımı.

Alt modüller:
    optics       kapalı-form türbülans bağıntıları (doğrulama referansı)
    screen       Kolmogorov faz ekranı (FFT + subharmonik)
    fields       tilt / PSF / sintilasyon alanları
    simulator    kare-tutarlı bozulma simülatörü (TurbParams, TurbulenceSimulator)
    hooks        mevcut ön-uca minimal enjeksiyon (türbülans kapalıyken birebir baseline)
    covariance   izoplanatik yama gruplaması ve Sigma_yama üreteci
    factors      GTSAM yama-gruplu ölçüm faktörü (A3/A4 kolları)
    metrics      ATE/RPE/NEES ve ön-uç metrikleri
    stats        eşleştirilmiş istatistik üçlüsü (bootstrap + Wilcoxon + işaret)
"""
from .simulator import TurbParams, TurbulenceSimulator  # noqa: F401

__all__ = ["TurbParams", "TurbulenceSimulator"]
