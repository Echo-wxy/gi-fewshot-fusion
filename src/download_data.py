"""Resumable data downloader. Usage: python src/download_data.py [--only endo,stad,crc,tcga]
Root folder comes from env var GIFS_ROOT (default ~/gifs_data).
Files you already downloaded with a browser are detected and verified, not downloaded again."""
import argparse, hashlib, os, sys, tarfile, time, zipfile
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("the requests package is missing: pip install -r requirements.txt")
try:  # verify HTTPS with the operating-system certificate store instead of certifi
    import truststore
    truststore.inject_into_ssl()
except ImportError:
    pass

SSL_HELP = ("    HTTPS certificate verification failed (often caused by a VPN or proxy).\n"
            "    Retry without the VPN or proxy, or download the file in a browser from the URL in the summary.")

ROOT = Path(os.environ.get("GIFS_ROOT", str(Path.home() / "gifs_data")))
DATA = ROOT / "data"
Z = "https://zenodo.org/api/records/2530835/files/{}/content"
GH = "https://github.com/cBioPortal/datahub/raw/e1b75cd25bf7436d16719b7d76117ad8f3f66e54/public/{}/{}"

# group, file name, url, sub-folder, exact size (bytes), md5, post-processing
PAGES = {"endo": "https://datasets.simula.no/hyper-kvasir/  (file hyper-kvasir-labeled-images.zip in the download table)",
         "stad": "https://zenodo.org/records/2530835", "crc": "https://zenodo.org/records/2530835",
         "tcga": "https://github.com/cBioPortal/datahub/tree/e1b75cd25bf7436d16719b7d76117ad8f3f66e54/public  (open the cohort folder, click a file name, then Download raw file)"}

ITEMS = [
    ("endo", "hyper-kvasir-labeled-images.zip",
     "https://datasets.simula.no/downloads/hyper-kvasir/hyper-kvasir-labeled-images.zip",
     "endoscopy", None, None, "unzip"),
    ("stad", "STAD_TRAIN_MSS.zip", Z.format("STAD_TRAIN_MSS.zip"), "wsi_tiles/STAD", 5848305951, "273295ed37ebccc3a574ab9106145e32", None),
    ("stad", "STAD_TRAIN_MSIMUT.zip", Z.format("STAD_TRAIN_MSIMUT.zip"), "wsi_tiles/STAD", 5841008485, "c785506ccdf49d4d08b65c5e7cdd28ae", None),
    ("stad", "STAD_TEST_MSS.zip", Z.format("STAD_TEST_MSS.zip"), "wsi_tiles/STAD", 10484110033, "8f324496c6feffe5353f493ae5e8c9ff", None),
    ("stad", "STAD_TEST_MSIMUT.zip", Z.format("STAD_TEST_MSIMUT.zip"), "wsi_tiles/STAD", 3226442356, "77941b8914a2fac469ce4843c75ee19a", None),
    ("crc", "CRC_DX_TRAIN_MSS.zip", Z.format("CRC_DX_TRAIN_MSS.zip"), "wsi_tiles/CRC", 5259435967, "0dbb048590925920c383e423fa282f77", None),
    ("crc", "CRC_DX_TRAIN_MSIMUT.zip", Z.format("CRC_DX_TRAIN_MSIMUT.zip"), "wsi_tiles/CRC", 5256478484, "d98f4ee32f5923a09779d9aabd773e37", None),
    ("crc", "CRC_DX_TEST_MSS.zip", Z.format("CRC_DX_TEST_MSS.zip"), "wsi_tiles/CRC", 7890181037, "3407a53dddec68e67f49b863a3520e88", None),
    ("crc", "CRC_DX_TEST_MSIMUT.zip", Z.format("CRC_DX_TEST_MSIMUT.zip"), "wsi_tiles/CRC", 3251312941, "006962fd57ffa71db43f3d768a63d9ec", None),
    # cBioPortal datahub on GitHub (git-LFS), pinned to commit e1b75cd (2026-10-30); sha256 from the LFS pointers
    ("tcga", "data_clinical_patient.txt", GH.format("stad_tcga_pan_can_atlas_2018", "data_clinical_patient.txt"), "tcga/stad_tcga_pan_can_atlas_2018", 134386, "sha256:e92acfd65317e19e8682460999f984f6b57d340dc38ed184bce3b6fa8411a4d5", None),
    ("tcga", "data_clinical_sample.txt", GH.format("stad_tcga_pan_can_atlas_2018", "data_clinical_sample.txt"), "tcga/stad_tcga_pan_can_atlas_2018", 93194, "sha256:e082f706bbbb5ff1104389b93cdc852f9dd31019b0bcf4afcbf7d7654cb2ace1", None),
    ("tcga", "data_mrna_seq_v2_rsem.txt", GH.format("stad_tcga_pan_can_atlas_2018", "data_mrna_seq_v2_rsem.txt"), "tcga/stad_tcga_pan_can_atlas_2018", 134512447, "sha256:cf741a33df0c0f2783e7d84de822bb2104b7fc80bd8423d64e9deb8db9a01217", None),
    ("tcga", "data_cna.txt", GH.format("stad_tcga_pan_can_atlas_2018", "data_cna.txt"), "tcga/stad_tcga_pan_can_atlas_2018", 24627999, "sha256:4780337064e69c0f1d68a1338fb5560163479674d4328d3339e8601aa93e742d", None),
    ("tcga", "data_mutations.txt", GH.format("stad_tcga_pan_can_atlas_2018", "data_mutations.txt"), "tcga/stad_tcga_pan_can_atlas_2018", 237607789, "sha256:84da71c9c20598ebc5e62a8ee078466d780abe31d5d2415918058a14c77bf808", None),
    ("tcga", "data_clinical_patient.txt", GH.format("coadread_tcga_pan_can_atlas_2018", "data_clinical_patient.txt"), "tcga/coadread_tcga_pan_can_atlas_2018", 175194, "sha256:3c5738f41b9047313af8f809c244bc66419991379b0c1978816dfc56c184c313", None),
    ("tcga", "data_clinical_sample.txt", GH.format("coadread_tcga_pan_can_atlas_2018", "data_clinical_sample.txt"), "tcga/coadread_tcga_pan_can_atlas_2018", 105858, "sha256:b9cbcc7cbc7471b44a3767aac48d4d0625001e3c07f8776c50f77364341f7817", None),
    ("tcga", "data_mrna_seq_v2_rsem.txt", GH.format("coadread_tcga_pan_can_atlas_2018", "data_mrna_seq_v2_rsem.txt"), "tcga/coadread_tcga_pan_can_atlas_2018", 100912131, "sha256:3e297f6ff45f7e21b1401b8c49b657ef2f1037728b9edf1e65a13829f680dcc4", None),
    ("tcga", "data_cna.txt", GH.format("coadread_tcga_pan_can_atlas_2018", "data_cna.txt"), "tcga/coadread_tcga_pan_can_atlas_2018", 32664387, "sha256:ae5e25e0017450326e5ff3f11aceb388b1d8ace056a3b86bc5c0ae8e105cf2c9", None),
    ("tcga", "data_mutations.txt", GH.format("coadread_tcga_pan_can_atlas_2018", "data_mutations.txt"), "tcga/coadread_tcga_pan_can_atlas_2018", 342402448, "sha256:f3e27a06e4497767bdae99e993a4f627493b6cce1f7842658b102bd6dcc7091a", None),]


def human(n):
    return f"{n / 1e9:.2f} GB" if n >= 1e8 else f"{n / 1e6:.1f} MB"


def md5_of(path, chunk=8 << 20, algo="md5"):
    h = hashlib.new(algo)
    done, total, t0 = 0, path.stat().st_size, time.time()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
            done += len(b)
            if time.time() - t0 > 5:
                print(f"\r    verifying {done * 100 // total}%", end="", flush=True)
                t0 = time.time()
    print("\r", end="")
    return h.hexdigest()


def probe(url):
    """Ask for byte 0 only. Returns (total_size, supports_range);
    (status_code, None) on an HTTP error; (None, False) if unreachable."""
    try:
        h = {"Accept-Encoding": "identity", "Range": "bytes=0-0"}
        with requests.get(url, headers=h, stream=True, timeout=30, allow_redirects=True) as r:
            if r.status_code >= 400 and r.status_code != 416:
                return r.status_code, None
            if r.status_code == 206 and "/" in r.headers.get("Content-Range", ""):
                return int(r.headers["Content-Range"].rsplit("/", 1)[1]), True
            return int(r.headers.get("Content-Length", 0)) or None, False
    except requests.exceptions.SSLError:
        return "ssl", None
    except Exception:
        return None, False


def _fetch_range(url, seg_file, start, end, stop, errors, retries=60):
    """Download bytes [start, end] into seg_file, resuming from its current size."""
    need = end - start + 1
    ssl_fail = 0
    for _ in range(retries):
        if stop.is_set():
            return
        have = seg_file.stat().st_size if seg_file.exists() else 0
        if have > need:
            seg_file.unlink()
            have = 0
        if have == need:
            return
        try:
            h = {"Accept-Encoding": "identity", "Range": f"bytes={start + have}-{end}"}
            with requests.get(url, headers=h, stream=True, timeout=60) as r:
                if r.status_code == 429 or r.status_code >= 500:
                    time.sleep(30)
                    continue
                if r.status_code != 206:
                    errors.append(f"server returned {r.status_code}")
                    stop.set()
                    return
                with open(seg_file, "ab") as f:
                    for chunk in r.iter_content(1 << 20):
                        if stop.is_set():
                            return
                        f.write(chunk)
        except requests.exceptions.SSLError:
            ssl_fail += 1
            if ssl_fail >= 3:
                errors.append("HTTPS certificate verification failed")
                stop.set()
                return
            time.sleep(5)
        except Exception:
            time.sleep(5)
    if not (seg_file.exists() and seg_file.stat().st_size == need):
        errors.append("retries exhausted")


def download_parallel(url, dst, size, n_conn=8, chunk=64 << 20):
    """Zenodo serves ~1 MB/s per connection, so fetch 64 MB pieces with 8 connections at once.
    Each piece is its own file (dst.pNNNN); a rerun skips finished pieces and resumes partial ones."""
    import queue, threading
    pieces = [(i, i * chunk, min(size, (i + 1) * chunk) - 1) for i in range(-(-size // chunk))]
    files = {i: dst.with_name(f"{dst.name}.p{i:04d}") for i, _, _ in pieces}
    todo = queue.Queue()
    for pc in pieces:
        todo.put(pc)
    stop, errors = threading.Event(), []

    def worker():
        while not stop.is_set():
            try:
                i, a_, b_ = todo.get_nowait()
            except queue.Empty:
                return
            _fetch_range(url, files[i], a_, b_, stop, errors)

    threads = [threading.Thread(target=worker, daemon=True) for _ in range(n_conn)]
    for t in threads:
        t.start()
    cur = lambda: sum(f.stat().st_size for f in files.values() if f.exists())
    t0, got0 = time.time(), cur()
    try:
        while any(t.is_alive() for t in threads):
            time.sleep(2)
            got = cur()
            speed = (got - got0) / max(time.time() - t0, 1e-6)
            print(f"\r    {got * 100 / size:5.1f}%  {human(got)} / {human(size)}  {speed / 1e6:.1f} MB/s ({n_conn} parallel connections)   ",
                  end="", flush=True)
    except KeyboardInterrupt:
        stop.set()
        raise
    print()
    if errors or cur() != size:
        print(f"    incomplete ({errors[0] if errors else 'some chunks interrupted'}); run again to resume")
        return False
    print("    merging chunks ...")
    with open(dst, "wb") as out:
        for i, _, _ in pieces:
            with open(files[i], "rb") as src:
                while True:
                    blk = src.read(16 << 20)
                    if not blk:
                        break
                    out.write(blk)
    for f in files.values():
        f.unlink()
    return True


def download(url, dst, size=None, retries=30):
    probed, ranges = probe(url)
    if probed == "ssl":
        print(SSL_HELP)
        return False
    if ranges is None and probed and 400 <= probed < 500 and probed != 429:
        print(f"    server returned {probed}; this URL cannot be downloaded directly")
        return False
    total = size or (probed if ranges is not None else None)
    if ranges and total and total > 200_000_000:
        return download_parallel(url, dst, total)
    part = dst.with_suffix(dst.suffix + ".part")
    ssl_count = 0
    for attempt in range(1, retries + 1):
        have = part.stat().st_size if part.exists() else 0
        if size and have == size:
            break
        if size and have > size:          # corrupt partial file: start over
            part.unlink()
            have = 0
        headers = {"Accept-Encoding": "identity"}
        if have:
            headers["Range"] = f"bytes={have}-"
        try:
            with requests.get(url, headers=headers, stream=True, timeout=60) as r:
                if r.status_code == 416:      # nothing left to fetch
                    break
                r.raise_for_status()
                if have and r.status_code != 206:   # server ignored Range: restart
                    have = 0
                total = size or (have + int(r.headers.get("Content-Length", 0)))
                mode = "ab" if have else "wb"
                t0, last, got0 = time.time(), time.time(), have
                with open(part, mode) as f:
                    for chunk in r.iter_content(4 << 20):
                        f.write(chunk)
                        have += len(chunk)
                        if time.time() - last > 2:
                            speed = (have - got0) / max(time.time() - t0, 1e-6)
                            pct = f"{have * 100 / total:5.1f}%" if total else ""
                            print(f"\r    {pct} {human(have)} / {human(total) if total else '?'}  {speed / 1e6:.1f} MB/s   ",
                                  end="", flush=True)
                            last = time.time()
            print()
            if size is None or have == size:
                break
        except requests.HTTPError as e:
            code = e.response.status_code if e.response is not None else 0
            if 400 <= code < 500 and code != 429:
                print(f"\n    server returned {code}; this URL cannot be downloaded directly")
                return False
            print(f"\n    server error {code} (attempt {attempt}), retrying in 10 s")
            time.sleep(10)
        except requests.exceptions.SSLError:
            ssl_count += 1
            if ssl_count >= 3:
                print("\n" + SSL_HELP)
                return False
            print(f"\n    HTTPS certificate verification failed (attempt {ssl_count}), retrying in 10 s")
            time.sleep(10)
        except Exception as e:  # network hiccup: retry and resume
            print(f"\n    connection interrupted (attempt {attempt}): {type(e).__name__}, resuming in 10 s")
            time.sleep(10)
    else:
        return False
    part.replace(dst)
    return True


def verify(path, size, md5):
    ok_flag = path.with_suffix(path.suffix + ".ok")
    if ok_flag.exists():
        # the flag stores the size the file had when it passed; any later change forces a full re-check
        if ok_flag.read_text().strip() == str(path.stat().st_size):
            return True
        ok_flag.unlink()
    if size and path.stat().st_size != size:
        print(f"    size mismatch: {path.stat().st_size} != {size}")
        return False
    if md5:
        algo, want = md5.split(":", 1) if ":" in md5 else ("md5", md5)
        got = md5_of(path, algo=algo)
        if got != want:
            print(f"    {algo} mismatch: {got}")
            return False
    else:  # no published checksum: check the archive opens
        try:
            if path.suffix == ".zip":
                zipfile.ZipFile(path).namelist()
            else:
                with tarfile.open(path) as t:
                    t.getmembers()
        except Exception as e:
            print(f"    cannot open archive: {e}")
            return False
    ok_flag.write_text(str(path.stat().st_size))
    return True


def post(path, how, out_dir):
    flag = out_dir / f".{path.name}.extracted"
    if flag.exists() or how is None:
        return
    print(f"    extracting to {out_dir} ...")
    if how == "unzip":
        zipfile.ZipFile(path).extractall(out_dir)
    else:
        with tarfile.open(path) as t:
            t.extractall(out_dir)
    flag.write_text("ok")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="endo,stad,crc,tcga")
    groups = set(ap.parse_args().only.split(","))
    print(f"data root: {DATA}")
    failed = []
    for g, name, url, sub, size, md5, how in ITEMS:
        if g not in groups:
            continue
        out = DATA / sub
        out.mkdir(parents=True, exist_ok=True)
        dst = out / name
        print(f"[{g}] {sub}/{name}" + (f"（{human(size)}）" if size else ""))
        if dst.exists() and verify(dst, size, md5):
            print("    exists and verified, skipping")
        else:
            if dst.exists():
                dst.rename(dst.with_suffix(dst.suffix + ".bad"))
                print("    existing file failed verification, renamed to .bad, downloading again")
            if not download(url, dst, size) or not verify(dst, size, md5):
                failed.append((name, PAGES[g], out))
                continue
            print("    downloaded and verified")
        post(dst, how, out)
    print("\n================ download summary ================")
    if not failed:
        print("all data ready.")
    else:
        print("The following files failed. Run the script again to resume,")
        print("or download them in a browser, keep the file names, place them in the folders below and run the script again:")
        for name, url, out in failed:
            print(f"  {name}\n    URL: {url}\n    save to: {out}")
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
