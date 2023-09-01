**安装**:
- [Biopython](https://biopython.org/wiki/Download)
```
pip install biopython --upgrade
```

**使用**:
```
python main.py
```

**样例**:
在comp_before.pdb 中检索`OD1 ASP R 121`, 然后找到对应的坐标，再用坐标在comp.pdb中找到 `OD1 ASP    79`
输入` OD1 ASP R 121` 输出 `OD1 ASP    79`