# -*- coding: utf-8 -*-
r"""
cbct_cut_1.py —— 按目录名配对执行 CT 切割

相对原始 cbct_cut_1.py 的改动（只改文件发现/配对逻辑，切割与可视化算法完全保留）：
  1. 原脚本用 glob 在 nii_dir / txt_dir 平铺查找文件，实际数据是嵌套结构：
         <nii_root>\\<前缀>\\*.nii.gz
         <txt_root>\\<前缀>rabbit\\*.txt
     现在改为：取目录名前三位数字作为匹配键，两侧键相同的目录配成一组。
  2. 组内按文件名排序逐个配对（nii 数量与 txt 数量必须相等，否则跳过该组并报错）。
  3. 输出目录：ct_cut\\<前缀>rabbit\\<txt文件名>\\<序号>.nii.gz
                ct_cut\\<前缀>rabbit\\visualizations\\<txt文件名>_visualization.png
  4. 默认"跳过已存在的输出文件"：逐对检查，缺什么补什么；整对齐全时直接跳过（不加载 NIfTI，
     大幅省时）。加 --force 可全部重算覆盖。
  5. 可用 --out 指定输出根目录（默认 ct_cut），例如：
         python cbct_cut_1.py --out <输出根目录>
  6. 【尺寸溢出修复】原写法 x_max = int(x + 宽/2) 与 x_min = int(x - 宽/2) 各自独立取整，
     当 x 的小数部分恰好卡在浮点边界时（例如 243.99999999999997），两者之差会比边长多 1，
     切出 161 格的 Y 方向立方体。现改为先取整边长，再令 max = min + 边长，
     保证 X 永远不超过 120、Y 永远不超过 160（贴边裁剪只会更小）。
  7. 可用 --fix-overflow 只扫描并重切已存在的溢出立方体（就地替换），不必整体重跑：
         python cbct_cut_1.py --out <输出根目录> --fix-overflow
  8. 数据集目录布局与配对键均可配置（默认沿用旧数据集）：
         --nii-root <NIfTI 根目录>   --txt-root <TXT 根目录>
         --match prefix3  ：目录名前三位数字相同即配对（旧数据集：ct_nii\809 ↔ set\809rabbit）
         --match dirname  ：递归查找，目录名完全相同才配对（新数据集：nii\rabbit_ct_nii\621rabbit1
                            ↔ set_all\621rabbit1）
     输出子目录名统一取 **TXT 目录名**，例如 809rabbit / 621rabbit1，两种模式都自然成立。
"""

import os
import sys
import glob
import argparse
import gc
import traceback

import numpy as np

try:  # 避免 Windows 控制台 GBK 编码导致的打印异常
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# ==================== 路径与参数 ====================
# Defaults only; every one of them is overridden by the command-line options below.
NII_ROOT = "ct_nii"      # NIfTI root (holding the <prefix> subdirectories)
TXT_ROOT = "set"         # coordinate-record root (holding the <prefix>rabbit subdirectories)
OUT_ROOT = "ct_cut"      # output root

PLANE_SIZE_PX = (640, 640)      # 平面像素尺寸
PLANE_SIZE_CM = (16.0, 16.0)    # 平面对应物理尺寸 (cm)
CUBE_SIZE_CM = (3, 4)           # 立方体尺寸 (cm)
Z_RANGE = (160, 300)            # Z 方向范围
PREFIX_LEN = 3                  # 目录名前三位数字作为配对键


# ==================== 输入解析 ====================
def load_and_prepare(nii_path):
    """加载NIfTI文件并保持原始方向"""
    img = nib.load(nii_path)
    return img, img.get_fdata()


def read_input_file(txt_path):
    """从文本文件读取参考点和目标点坐标（与原始脚本一致）"""
    targets = []
    ref_point_px = None
    ref_point_cm = None

    base_name = os.path.splitext(os.path.basename(txt_path))[0]

    with open(txt_path, 'r', encoding='utf-8', errors='replace') as f:
        first_line = f.readline().strip()
        if first_line and not first_line.startswith('#'):
            try:
                parts = first_line.split('#')[0].split()
                if len(parts) >= 2:
                    ref_point_px = (float(parts[0]), float(parts[1]))
                if len(parts) >= 4:
                    ref_point_cm = (float(parts[2]), float(parts[3]))
            except Exception:
                pass

        second_line = f.readline().strip()
        if not ref_point_cm and second_line and not second_line.startswith('#'):
            try:
                parts = second_line.split('#')[0].split()
                if len(parts) >= 2:
                    ref_point_cm = (float(parts[0]), float(parts[1]))
            except Exception:
                pass

        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            try:
                coords = [float(x.strip()) for x in line.replace(',', ' ').split() if x.strip()]
                if len(coords) >= 2:
                    targets.append((coords[0], coords[1]))
            except Exception as e:
                print(f"警告：跳过无法解析的行 '{line}': {str(e)}")

    if not ref_point_px:
        ref_point_px = (528, 294)
        print(f"警告：未找到有效的参考点像素坐标，使用默认值: {ref_point_px}")

    if not ref_point_cm:
        ref_point_cm = (15.0, 10.0)
        print(f"警告：未找到有效的参考点物理坐标，使用默认值: {ref_point_cm}")

    print(f"从文件读取: 参考点像素 = {ref_point_px}, 参考点物理 = {ref_point_cm}")
    print(f"目标点数量: {len(targets)}")

    return ref_point_px, ref_point_cm, targets, base_name


# ==================== 几何计算 ====================
def compute_boxes(target_list, ref_point_px, ref_point_cm,
                  plane_size_px, plane_size_cm, cube_size_cm,
                  z_range, z_max_limit):
    """按原始脚本的同一套公式计算每个目标点的像素中心与切割框（不读图像即可完成）"""
    px_per_cm_x = plane_size_px[0] / plane_size_cm[0]
    px_per_cm_y = plane_size_px[1] / plane_size_cm[1]

    # 边长先取整：X=120、Y=160。随后用 max = min + 边长 定上限，
    # 这样两个端点不会再各自取整而相差 1 格（浮点边界导致的 161 溢出）。
    cube_width_px = int(round(cube_size_cm[0] * px_per_cm_x))
    cube_height_px = int(round(cube_size_cm[1] * px_per_cm_y))

    boxes = []
    for target_cm in target_list:
        dx_cm = target_cm[0] - ref_point_cm[0]
        dy_cm = ref_point_cm[1] - target_cm[1]  # Y轴反转处理

        target_x = ref_point_px[0] + dx_cm * px_per_cm_x
        target_y = ref_point_px[1] + dy_cm * px_per_cm_y

        target_x_adjusted = target_x
        target_y_adjusted = target_y

        x_min = int(target_x_adjusted - cube_width_px / 2)
        x_max = x_min + cube_width_px
        y_min = int(target_y_adjusted - cube_height_px / 2)
        y_max = y_min + cube_height_px

        # 边界检查：只允许裁小，绝不允许超过 cube_width_px / cube_height_px
        x_min = max(0, x_min)
        x_max = min(plane_size_px[0], x_max)
        y_min = max(0, y_min)
        y_max = min(plane_size_px[1], y_max)
        if x_max - x_min > cube_width_px:
            x_max = x_min + cube_width_px
        if y_max - y_min > cube_height_px:
            y_max = y_min + cube_height_px
        z_min = max(0, z_range[0])
        z_max = min(z_max_limit, z_range[1])

        boxes.append({
            'bbox': (x_min, x_max, y_min, y_max),
            'z': (z_min, z_max),
            'target': target_cm,
            'center_px': (target_x_adjusted, target_y_adjusted),
            'original_center_px': (target_x, target_y),
        })
    return boxes


# ==================== 可视化（与原脚本一致） ====================
def visualize_results(data, results, z_range, plane_size_px, plane_size_cm,
                      ref_point_px, base_name, vis_dir):
    """完全匹配3D Slicer的可视化 - 使用txt文件名命名"""
    z_mid = (z_range[0] + z_range[1]) // 2
    slice_index = min(z_mid, data.shape[2] - 1)
    original_slice = data[:, :, slice_index]
    display_img = original_slice.T

    plt.figure(figsize=(12, 12))
    plt.imshow(display_img, cmap='gray',
               origin='upper',
               extent=[0, plane_size_px[0], plane_size_px[1], 0],
               vmin=-100, vmax=1000)

    plt.scatter(ref_point_px[0], ref_point_px[1],
                s=300, marker='*',
                color='red', edgecolor='white',
                zorder=10, alpha=0.9,
                label=f'参考点 ({ref_point_px[0]}, {ref_point_px[1]})')

    for i, res in enumerate(results):
        x_min, x_max, y_min, y_max = res['bbox']
        rect_width = x_max - x_min
        rect_height = y_max - y_min

        rect = Rectangle((x_min, y_min),
                         rect_width, rect_height,
                         edgecolor=plt.cm.tab20(i),
                         facecolor='none',
                         linewidth=2,
                         alpha=0.8)
        plt.gca().add_patch(rect)

        center_x = (x_min + x_max) / 2
        center_y = (y_min + y_max) / 2

        plt.text(center_x, center_y,
                 f"{i + 1}\n({res['target'][0]:.1f}, {res['target'][1]:.1f})",
                 color='white',
                 ha='center', va='center',
                 fontsize=10,
                 fontweight='bold',
                 bbox=dict(facecolor=plt.cm.tab20(i), alpha=0.7))

    scale_length_cm = 1.5
    scale_length_px = scale_length_cm * (plane_size_px[0] / plane_size_cm[0])

    plt.plot([plane_size_px[0] / 2 - scale_length_px / 2,
              plane_size_px[0] / 2 + scale_length_px / 2],
             [plane_size_px[1] - 50, plane_size_px[1] - 50],
             'w-', linewidth=4)
    plt.text(plane_size_px[0] / 2, plane_size_px[1] - 70,
             f"{scale_length_cm} cm",
             color='white', ha='center', va='top',
             fontsize=11)

    plt.text(plane_size_px[0] - 20, 30,
             f"ZE: {z_mid} (T:{z_range[0]}-{z_range[1]})",
             color='white', ha='right', va='top',
             fontsize=12)

    plt.text(10, plane_size_px[1] - 20,
             f"B: {base_name}",
             color='white', ha='left', va='bottom',
             fontsize=11)

    plt.title(f"{base_name} - Z层: {z_mid} (范围 {z_range[0]}-{z_range[1]})", fontsize=14)
    plt.legend(loc='lower right')
    plt.axis('off')

    vis_path = os.path.join(vis_dir, f"{base_name}_visualization.png")
    plt.tight_layout()
    plt.savefig(vis_path, dpi=150)
    plt.close()
    print(f"✅ 结果可视化已保存为 '{vis_path}'")
    return vis_path


# ==================== 配对发现 ====================
def prefix_key(name):
    """取目录名中最先出现的连续三位数字作为配对键（如 '809' / '809rabbit' -> '809'）"""
    digits = ""
    for ch in name:
        if ch.isdigit():
            digits += ch
            if len(digits) == PREFIX_LEN:
                return digits
        else:
            digits = ""
    return None


def discover_pairs(nii_root, txt_root, mode="prefix3"):
    """配对 NIfTI 目录与 TXT 目录。

    返回 (pairs, nii_only, txt_only)：
      pairs —— {组名: (nii目录, txt目录)}，组名 = TXT 目录名（同时用作输出子目录名）
      nii_only —— 只有 NIfTI 没有 TXT 的目录名列表
      txt_only —— 只有 TXT 没有 NIfTI 的目录名列表
    """
    if mode == "dirname":
        def scan(root, suffix):
            found, dup = {}, []
            for dirpath, _dirnames, filenames in os.walk(root):
                if any(f.endswith(suffix) for f in filenames):
                    key = os.path.basename(dirpath)
                    if key in found:
                        dup.append((key, dirpath))
                    else:
                        found[key] = dirpath
            return found, dup

        nii_found, nii_dup = scan(nii_root, ".nii.gz")
        txt_found, txt_dup = scan(txt_root, ".txt")
        for key, p in nii_dup + txt_dup:
            print(f"警告: 目录名 '{key}' 重复出现（{p}），只保留第一个")

        pairs = {k: (nii_found[k], txt_found[k])
                 for k in sorted(set(nii_found) & set(txt_found))}
        return pairs, sorted(set(nii_found) - set(txt_found)), sorted(set(txt_found) - set(nii_found))

    # ---- prefix3：目录名前三位数字相同即配对（旧数据集） ----
    nii_by_prefix = {}
    for name in sorted(os.listdir(nii_root)):
        path = os.path.join(nii_root, name)
        if os.path.isdir(path):
            key = prefix_key(name)
            if key and key not in nii_by_prefix:
                nii_by_prefix[key] = path

    pairs, txt_only, used = {}, [], set()
    for name in sorted(os.listdir(txt_root)):
        path = os.path.join(txt_root, name)
        if not os.path.isdir(path):
            continue
        key = prefix_key(name)
        if key and key in nii_by_prefix:
            pairs[name] = (nii_by_prefix[key], path)   # 组名 = TXT 目录名
            used.add(key)
        else:
            txt_only.append(name)

    nii_only = sorted(k for k in nii_by_prefix if k not in used)
    return pairs, nii_only, txt_only


# ==================== 单个文件对处理 ====================
def process_pair(nii_path, txt_path, group_name, force=False, out_root=None):
    """处理一个 (nii, txt) 对；返回统计字典。默认跳过已存在的输出文件。"""
    if out_root is None:
        out_root = OUT_ROOT
    txt_base = os.path.splitext(os.path.basename(txt_path))[0]
    group_dir = os.path.join(out_root, group_name)
    out_dir = os.path.join(group_dir, txt_base)
    vis_dir = os.path.join(group_dir, "visualizations")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(vis_dir, exist_ok=True)

    ref_point_px, ref_point_cm, target_list, base_name = read_input_file(txt_path)

    stats = {'cube_written': 0, 'cube_skipped': 0, 'vis_written': 0,
             'vis_skipped': 0, 'loaded': False}

    if not target_list:
        print(f"警告: {os.path.basename(txt_path)} 中没有找到有效的目标点，跳过")
        stats['empty'] = True
        return stats

    # 先不读图像，按公式算出每个目标点的框；z_max 上限先用 data.shape[2] 的占位，
    # 真正切割时再用实际体积形状裁剪（与原脚本 z_max = min(data.shape[2], z_range[1]) 等价）。
    boxes = compute_boxes(target_list, ref_point_px, ref_point_cm,
                          PLANE_SIZE_PX, PLANE_SIZE_CM, CUBE_SIZE_CM,
                          Z_RANGE, z_max_limit=Z_RANGE[1])

    cube_paths = [os.path.join(out_dir, f"{i + 1}.nii.gz") for i in range(len(boxes))]
    vis_path = os.path.join(vis_dir, f"{base_name}_visualization.png")

    pending = list(range(len(boxes))) if force else \
        [i for i, p in enumerate(cube_paths) if not os.path.exists(p)]
    need_vis = force or (not os.path.exists(vis_path))

    stats['cube_skipped'] = len(boxes) - len(pending)
    stats['vis_skipped'] = 0 if need_vis else 1

    if not pending and not need_vis:
        print(f"  [跳过] 输出已齐全（{len(boxes)} 个立方体 + 可视化）")
        return stats

    print(f"  待生成立方体: {len(pending)}/{len(boxes)}，可视化: {'需要' if need_vis else '已存在'}")

    img, data = load_and_prepare(nii_path)
    stats['loaded'] = True

    # 用真实体积形状修正 z 范围
    for b in boxes:
        z_min, z_max = b['z']
        b['z'] = (z_min, min(data.shape[2], z_max))
    z_range_eff = (max(0, Z_RANGE[0]), min(data.shape[2], Z_RANGE[1]))

    for i in pending:
        b = boxes[i]
        x_min, x_max, y_min, y_max = b['bbox']
        z_min, z_max = b['z']

        cube_data = data[x_min:x_max, y_min:y_max, z_min:z_max]
        save_cube(cube_data, img, cube_paths[i])
        stats['cube_written'] += 1
        print(f"  立方体 {i + 1} 已保存: {os.path.basename(cube_paths[i])} "
              f"(bbox={b['bbox']}, z={b['z']}, shape={cube_data.shape})")

    if need_vis:
        visualize_results(data, boxes, z_range_eff, PLANE_SIZE_PX, PLANE_SIZE_CM,
                          ref_point_px, base_name, vis_dir)
        stats['vis_written'] = 1

    del data, img
    gc.collect()
    return stats


# ==================== 尺寸上限与溢出修复 ====================
def cube_limits():
    """立方体允许的最大体素格数 (X, Y)"""
    return (int(round(CUBE_SIZE_CM[0] * PLANE_SIZE_PX[0] / PLANE_SIZE_CM[0])),
            int(round(CUBE_SIZE_CM[1] * PLANE_SIZE_PX[1] / PLANE_SIZE_CM[1])))


def save_cube(cube_data, img, path):
    """保存立方体，并强制校验不超出 X/Y 上限"""
    lim_x, lim_y = cube_limits()
    if cube_data.shape[0] > lim_x or cube_data.shape[1] > lim_y:
        raise RuntimeError(f"立方体尺寸溢出: {cube_data.shape} 超过上限 ({lim_x}, {lim_y})")
    nib.save(nib.Nifti1Image(cube_data, img.affine, img.header), path)


def scan_overflow(out_root):
    """扫描输出目录，返回 {(组名, txt基名): [溢出的序号...]}；组名 = 输出的一级子目录名"""
    lim_x, lim_y = cube_limits()
    bad = {}
    for p in sorted(glob.glob(os.path.join(out_root, "*", "*", "*.nii.gz"))):
        rel = os.path.relpath(p, out_root)
        parts = rel.split(os.sep)
        if len(parts) != 3:
            continue
        group, base, fname = parts
        try:
            sh = nib.load(p).shape
        except Exception as e:
            print(f"警告: 无法读取 {rel}: {e}")
            continue
        if sh[0] > lim_x or sh[1] > lim_y:
            bad.setdefault((group, base), []).append(int(fname.split('.')[0]))
    return bad


def regenerate_cubes(nii_path, txt_path, group_name, indices, out_root, refresh_vis=True):
    """按修正后的几何重新切指定序号的立方体并就地覆盖，可选重画可视化"""
    txt_base = os.path.splitext(os.path.basename(txt_path))[0]
    group_dir = os.path.join(out_root, group_name)
    out_dir = os.path.join(group_dir, txt_base)
    vis_dir = os.path.join(group_dir, "visualizations")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(vis_dir, exist_ok=True)

    ref_point_px, ref_point_cm, targets, base_name = read_input_file(txt_path)
    if not targets:
        print(f"  警告: {os.path.basename(txt_path)} 无有效目标点，跳过")
        return 0

    boxes = compute_boxes(targets, ref_point_px, ref_point_cm,
                          PLANE_SIZE_PX, PLANE_SIZE_CM, CUBE_SIZE_CM,
                          Z_RANGE, z_max_limit=Z_RANGE[1])

    img, data = load_and_prepare(nii_path)
    z_range_eff = (max(0, Z_RANGE[0]), min(data.shape[2], Z_RANGE[1]))

    n = 0
    for i in sorted(indices):
        if i < 1 or i > len(boxes):
            print(f"  警告: 序号 {i} 超出目标点范围，跳过")
            continue
        b = boxes[i - 1]
        x_min, x_max, y_min, y_max = b['bbox']
        z_min, z_max = b['z']
        z_max = min(data.shape[2], z_max)

        cube_data = data[x_min:x_max, y_min:y_max, z_min:z_max]
        dst = os.path.join(out_dir, f"{i}.nii.gz")
        save_cube(cube_data, img, dst)
        print(f"  已替换 {i}.nii.gz  bbox={b['bbox']}  z=({z_min}, {z_max})  shape={cube_data.shape}")
        n += 1

    if refresh_vis:
        visualize_results(data, boxes, z_range_eff, PLANE_SIZE_PX, PLANE_SIZE_CM,
                          ref_point_px, base_name, vis_dir)

    del data, img
    gc.collect()
    return n


def fix_overflow(out_root, refresh_vis=True, nii_root=None, txt_root=None, mode="prefix3"):
    """只重切已溢出的立方体并就地替换（同时重画受影响的配对可视化）"""
    nii_root = nii_root or NII_ROOT
    txt_root = txt_root or TXT_ROOT
    lim_x, lim_y = cube_limits()
    print(f"===== 溢出修复模式 =====")
    print(f"输出根目录: {out_root}")
    print(f"尺寸上限: X <= {lim_x}, Y <= {lim_y}（Z 不限制）")

    bad = scan_overflow(out_root)
    if not bad:
        print("未发现溢出立方体，无需修复。")
        return

    total_files = sum(len(v) for v in bad.values())
    print(f"发现溢出立方体 {total_files} 个，分布在 {len(bad)} 个文件对中\n")

    pairs, _nii_only, _txt_only = discover_pairs(nii_root, txt_root, mode)
    fixed = 0
    vis_n = 0
    failed = []

    for (group, base), idxs in sorted(bad.items()):
        print(f"\n{'-' * 66}")
        print(f"[修复] {group} / {base}: 溢出序号 {sorted(idxs)}")

        if group not in pairs:
            failed.append((group, base, "找不到对应的 NIfTI/TXT 目录"))
            print("  ❌ 找不到对应的源目录，跳过")
            continue
        nii_dir, txt_dir = pairs[group]

        nii_files = sorted(glob.glob(os.path.join(nii_dir, "*.nii.gz")))
        txt_files = sorted(glob.glob(os.path.join(txt_dir, "*.txt")))
        txt_path = os.path.join(txt_dir, base + ".txt")
        if txt_path not in txt_files:
            failed.append((group, base, "找不到对应 TXT"))
            print("  ❌ 找不到对应 TXT，跳过")
            continue
        i = txt_files.index(txt_path)
        if i >= len(nii_files):
            failed.append((group, base, "NIfTI 数量不足"))
            print("  ❌ NIfTI 数量不足，跳过")
            continue
        nii_path = nii_files[i]
        print(f"  源 NIfTI: {os.path.basename(nii_path)}")

        try:
            n = regenerate_cubes(nii_path, txt_path, group, idxs, out_root,
                                 refresh_vis=refresh_vis)
            fixed += n
            if refresh_vis:
                vis_n += 1
        except Exception as e:
            failed.append((group, base, str(e)))
            print(f"  ❌ 修复失败: {e}")
            traceback.print_exc()

    print(f"\n{'=' * 66}")
    print(f"修复完成: 重切并替换 {fixed} 个立方体，重画 {vis_n} 张可视化")
    print(f"失败条目: {len(failed)}")
    for group, base, msg in failed:
        print(f"  ! {group}/{base}: {msg}")


# ==================== 主流程 ====================
def main():
    parser = argparse.ArgumentParser(description="按目录名配对执行 CT 切割")
    parser.add_argument("--force", action="store_true",
                        help="忽略已存在的输出，全部重算覆盖")
    parser.add_argument("--prefix", default=None,
                        help="只处理指定组，逗号分隔（前缀匹配，例如 --prefix 812,621）")
    parser.add_argument("--out", default=OUT_ROOT,
                        help=f"输出根目录（默认 {OUT_ROOT}）")
    parser.add_argument("--nii-root", default=NII_ROOT,
                        help=f"NIfTI 根目录（默认 {NII_ROOT}）")
    parser.add_argument("--txt-root", default=TXT_ROOT,
                        help=f"TXT 根目录（默认 {TXT_ROOT}）")
    parser.add_argument("--match", choices=("prefix3", "dirname"), default="prefix3",
                        help="配对方式：prefix3=目录名前三位数字相同；dirname=目录名完全相同（递归查找）")
    parser.add_argument("--fix-overflow", action="store_true",
                        help="只扫描并重切已存在的溢出立方体（X>120 或 Y>160），就地替换")
    parser.add_argument("--no-vis", action="store_true",
                        help="配合 --fix-overflow：不重画可视化图片")
    args = parser.parse_args()

    out_root = os.path.abspath(args.out)
    nii_root = os.path.abspath(args.nii_root)
    txt_root = os.path.abspath(args.txt_root)
    only = set(p.strip() for p in args.prefix.split(',')) if args.prefix else None

    if args.fix_overflow:
        fix_overflow(out_root, refresh_vis=not args.no_vis,
                     nii_root=nii_root, txt_root=txt_root, mode=args.match)
        return

    print("===== 医学图像批量处理系统（按目录名配对版） =====", flush=True)
    print(f"NIfTI根目录: {nii_root}")
    print(f"TXT根目录:   {txt_root}")
    print(f"输出根目录:  {out_root}")
    print(f"配对方式:    {args.match}")
    print(f"模式:        {'全部重算覆盖 (--force)' if args.force else '跳过已存在的输出文件'}")

    pairs, only_nii, only_txt = discover_pairs(nii_root, txt_root, args.match)

    print(f"\n成功配对 {len(pairs)} 组")
    if only_nii:
        print(f"警告: 以下 NIfTI 目录没有对应的 TXT 目录，已忽略: {', '.join(only_nii)}")
    if only_txt:
        print(f"警告: 以下 TXT 目录没有对应的 NIfTI 目录，已忽略: {', '.join(only_txt)}")

    groups = sorted(pairs)
    if only:
        groups = [g for g in groups if any(g == o or g.startswith(o) for o in only)]
    print(f"本次处理 {len(groups)} 组: {', '.join(groups)}\n")

    total = {'pairs': 0, 'cube_written': 0, 'cube_skipped': 0,
             'vis_written': 0, 'vis_skipped': 0, 'failed': 0, 'empty': 0}
    summary = []

    for group in groups:
        nii_dir, txt_dir = pairs[group]

        nii_files = sorted(glob.glob(os.path.join(nii_dir, "*.nii.gz")))
        txt_files = sorted(glob.glob(os.path.join(txt_dir, "*.txt")))

        print(f"\n{'=' * 70}")
        print(f"组 {group}: NIfTI {len(nii_files)} 个, TXT {len(txt_files)} 个")

        if len(nii_files) != len(txt_files):
            print(f"❌ 错误: 组 {group} 的 NIfTI({len(nii_files)}) 与 TXT({len(txt_files)}) "
                  f"数量不匹配，跳过该组")
            summary.append((group, f"数量不匹配 NII={len(nii_files)} TXT={len(txt_files)}"))
            total['failed'] += 1
            continue
        if not nii_files:
            print(f"⚠️ 组 {group} 目录为空，跳过")
            summary.append((group, "空目录"))
            continue

        p_cube_w = p_cube_s = p_vis_w = p_vis_s = 0

        for i in range(len(nii_files)):
            nii_path = nii_files[i]
            txt_path = txt_files[i]

            print(f"\n{'-' * 70}")
            print(f"[{group}] 文件对 {i + 1}/{len(nii_files)}:")
            print(f"  NIfTI: {os.path.basename(nii_path)}")
            print(f"  TXT:   {os.path.basename(txt_path)}")

            total['pairs'] += 1
            try:
                st = process_pair(nii_path, txt_path, group, force=args.force,
                                  out_root=out_root)
                p_cube_w += st['cube_written']
                p_cube_s += st['cube_skipped']
                p_vis_w += st['vis_written']
                p_vis_s += st['vis_skipped']
                if st.get('empty'):
                    total['empty'] += 1
            except Exception as e:
                total['failed'] += 1
                print(f"❌ 处理失败: {str(e)}")
                traceback.print_exc()

        total['cube_written'] += p_cube_w
        total['cube_skipped'] += p_cube_s
        total['vis_written'] += p_vis_w
        total['vis_skipped'] += p_vis_s

        summary.append((group, f"立方体 新写{p_cube_w}/跳过{p_cube_s}, "
                               f"可视化 新写{p_vis_w}/跳过{p_vis_s}"))
        print(f"\n[{group}] 小计: 立方体 新写 {p_cube_w} / 跳过 {p_cube_s}; "
              f"可视化 新写 {p_vis_w} / 跳过 {p_vis_s}")

    print(f"\n{'=' * 70}")
    print("===== 总结 =====")
    for group, line in summary:
        print(f"  {group}: {line}")
    print(f"\n总文件对: {total['pairs']}")
    print(f"立方体: 新写 {total['cube_written']}, 跳过已存在 {total['cube_skipped']}")
    print(f"可视化: 新写 {total['vis_written']}, 跳过已存在 {total['vis_skipped']}")
    print(f"失败文件对: {total['failed']}, 无有效目标点: {total['empty']}")
    print("\n🎉 所有文件处理完成！")


if __name__ == "__main__":
    main()
