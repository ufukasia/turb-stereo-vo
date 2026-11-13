"""Her 10 s: surucu işçilerinin açık be.log dosyası (iş kimliği) ve birikimli CPU süresi (utime+stime)."""
import os, time, sys, glob
HZ = os.sysconf("SC_CLK_TCK")
out = open(sys.argv[1], "a")
while True:
    t = time.time()
    for p in glob.glob("/proc/[0-9]*"):
        try:
            cmd = open(p + "/cmdline", "rb").read()
            if b"spawn_main" not in cmd:
                continue
            job = None
            for fd in os.listdir(p + "/fd"):
                try:
                    l = os.readlink(f"{p}/fd/{fd}")
                except OSError:
                    continue
                if l.endswith("be.log") and "/r3090/fe/" in l:
                    job = l.split("/r3090/fe/")[1]
            if job is None:
                continue
            st = open(p + "/stat").read().rsplit(")", 1)[1].split()
            cpu = (int(st[11]) + int(st[12])) / HZ
            out.write(f"{t:.0f} {p[6:]} {cpu:.1f} {job}\n")
        except (OSError, IndexError):
            pass
    out.flush()
    time.sleep(10)
