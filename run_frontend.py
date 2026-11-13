"""
Frontend girişi: config.yaml'ı okuyup özellik çıkarımı + stereo eşleştirme + temporal takibi
çalıştırır, backend için kare başına landmark dump'ı üretir. Tüm ayarlar config.yaml'dan.

    python run_frontend.py        # -> outputs/landmarks_<seq>.pkl
"""
import os
import config
from stereo_front_end import StereoMatcher


def main(cfg=None):
    cfg = cfg or config.load()
    ds, f = cfg["dataset"], cfg["frontend"]
    os.makedirs(cfg["paths"]["outputs"], exist_ok=True)

    matcher = StereoMatcher(
        base_path=ds["base_path"],
        sequence=ds["sequence"],
        matcher_type=f["matcher"],
        filter_method=f["filter"],
        y_axis_threshold=f["y_axis_threshold"],
        extractor=f["extractor"],
        max_track_length=f["max_track_length"],
        max_keypoints=f["max_keypoints"],
        detection_threshold=f["detection_threshold"],
        stereo_extrinsics=f.get("stereo_extrinsics", "gomulu"),
    )
    config.apply_frontend(matcher, cfg)

    out = config.dump_path(cfg)
    matcher.process_sequence(max_frames=ds["max_frames"],
                             visualize=f["visualize"],
                             viz_interval=f["viz_interval"],
                             save_visualizations=f["save_visualizations"],
                             export_landmarks_path=out)
    # [TÜRBÜLANS] görüntü kalitesi özeti (PSNR/SSIM); türbülans kapalıyken dosya üretilmez
    if getattr(matcher, "turb", None) is not None:
        q = matcher.turb.quality_summary()
        if q is not None:
            import json
            with open(os.path.join(cfg["paths"]["outputs"], "turb_quality.json"), "w") as fq:
                json.dump(q, fq, indent=1)
    stats = matcher.compute_statistics()
    print(f"\n✓ Dump hazır: {out}")
    return stats


if __name__ == "__main__":
    main()
