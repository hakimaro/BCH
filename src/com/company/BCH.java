package com.company;

import java.util.ArrayList;

public class BCH {
    ArrayList<Integer> alpha;
    // Набор циклотомических классов
    ArrayList<ArrayList<Integer>> cycle = new ArrayList<>();
    // Порождающий многочлен g(x)
    boolean[] gx;
    // Примитивный многочлен
    ArrayList<boolean[]> reducingPolies = new ArrayList<>();
    boolean[] reducingPoly;
    // Мультипликативная группа
    ArrayList<boolean[]> mulGroup = new ArrayList<>();
    int m;
    int n;
    int d;

    /**
     * Инициализация БЧХ кода
     * @param m - степень
     * @param d -
     */
    public BCH(int m, int d) {
        this.m = m;
        this.n = (int) (Math.pow(2, m) - 1);
        this.d = d;
    }

    /**
     * Инициализация циклотомических классов
     */
    public void initCycle() {
        boolean[] used = new boolean[n + 1];
        cycle.add(new ArrayList<>());
        int i = findFirstUnused(used);
        while(findFirstUnused(used) != 0) {
            if (!used[i]) {
                used[i] = true;
                cycle.get(cycle.size() - 1).add(i);
                i = (i * 2) % n;
            } else {
                cycle.add(new ArrayList<>());
                i = findFirstUnused(used);
            }
        }
    }

    /**
     * Поиск элемента без цикла
     * @param used
     * @return
     */
    public int findFirstUnused(boolean[] used) {
        for(int i = 1; i <= n; i++) {
            if (!used[i]) return i;
        }
        return 0;
    }

    /**
     * Инициализация примитивного многочлена
     */
    public void initReducingPolinomial() {
        boolean[] dividend = new boolean[n + 1];
        dividend[n] = true; dividend[0] = true; // x^n - 1
        boolean[] divider = new boolean[m + 1];
        divider[m] = true; divider[0] = true;
        BinaryBruteForce(0, divider, dividend);
    }

    /**
     * Инициализация мультипликативной группы
     * @param reducingPoly
     */
    public void initMulGroup(boolean[] reducingPoly) {
        boolean[] mask = new boolean[getMaxDegree(reducingPoly)];
        System.arraycopy(reducingPoly, 0, mask, 0, getMaxDegree(reducingPoly));
        for(int i = 0; i < n; i++) {
            boolean[] alpha_i = new boolean[getMaxDegree(reducingPoly)];
            if (i < getMaxDegree(reducingPoly)) {
                alpha_i[i] = true;
            } else if (i == getMaxDegree(reducingPoly)) {
                System.arraycopy(mask, 0, alpha_i, 0, alpha_i.length);
            } else {
                for(int j = 0; j < mask.length; j++) {
                    if (mask[j]) {
                        alpha_i = sum(alpha_i, mulGroup.get(i - mask.length + j));
                    }
                }
            }
            mulGroup.add(alpha_i);
        }
    }

    /**
     * Вычисление порождающего многочлена
     */
    public void initGx() {
        boolean[] usedCycle = new boolean[cycle.size()];
        for(int alpha = 1; alpha < d; alpha++) {
            int cur_cycle = findCycle(alpha);
            if (usedCycle[cur_cycle]) continue;
            usedCycle[cur_cycle] = true;
            boolean[] f = new boolean[cycle.get(cur_cycle).size() + 1];
            f[f.length-1] = true;
            for(int degree = 0; degree < cycle.get(cur_cycle).size(); degree++) {
                boolean[] res = new boolean[1];
                this.alpha = new ArrayList<>();
                findAlpha(0, degree, 0, 0, cur_cycle);
                for(int i = 0; i < this.alpha.size(); i++) {
                    if (i == 0) res = mulGroup.get((this.alpha.get(i)) % mulGroup.size()).clone();
                    else res = sum(res, mulGroup.get((this.alpha.get(i)) % mulGroup.size())).clone();
                }
                f[degree] = res[0];
            }
            if (cur_cycle == 0) gx = f.clone();
            else gx = multiplying(gx, f).clone();
        }
    }

    /**
     * Вывод полинома
     * @param polynom
     */
    public static void printPolynom(boolean[] polynom) {
        if (getMaxDegree(polynom) == -1) System.out.print("0");
        for(int i = getMaxDegree(polynom); i >= 0; i--) {
            if (polynom[i]) System.out.print("x^" + i + ((i > 0) ? "+" : ""));
        }
        System.out.println();
    }

    /**
     * Поиск индекса циклотомического класса, в который входит альфа
     * @param alpha
     * @return
     */
    public int findCycle(int alpha) {
        for(int i = 0; i < this.cycle.size(); i++) {
            if (cycle.get(i).contains(alpha) ) return i;
        }
        return -1;
    }
    /**
     * Вывод циклотомических классов
     */
    public void printCycle() {
        for (ArrayList<Integer> cycles : cycle) {
            for (Integer num : cycles) {
                System.out.print(num + " ");
            }
            System.out.println();
        }
    }

    /**
     * Вывод мультипликативной группы
     */
    public void printMulGroup() {
        for(int i = 0; i < mulGroup.size(); i++) {
            System.out.print(i + ": ");
            printPolynom(mulGroup.get(i));
        }
    }

    /**
     * Деление многочлена по модулю 2
     * @param dividend -- делимое
     * @param divider -- делитель
     * @return remainder == 0
     */
    public boolean dividing(boolean[] dividend, boolean[] divider) {
        boolean[] remainder = dividend.clone();
        boolean[] tmp = new boolean[dividend.length]; // Промежуточные многочлены
        int r_remainder = getMaxDegree(remainder);
        int r_divider = getMaxDegree(divider);

        while (r_remainder >= r_divider) {
            if (r_divider + 1 >= 0)
                System.arraycopy(divider, 0, tmp, (r_remainder - r_divider), r_divider + 1);
            for(int i = r_remainder; i >= 0; i--) {
                remainder[i] ^= tmp[i];
            }
            r_remainder = getMaxDegree(remainder);
        }
        return getMaxDegree(remainder) == -1;
    }

    /** Умножение многочленов по модулю 2
     *
     * @param g1
     * @param g2
     * @return
     */
    public boolean[] multiplying(boolean[] g1, boolean[] g2) {
        boolean[] res = new boolean[getMaxDegree(g1) + getMaxDegree(g2) + 1];
        for(int i = 0; i <= getMaxDegree(g1); i++) {
            for(int j = 0; j <= getMaxDegree(g2); j++) {
                res[i + j] ^= (g1[i] & g2[j]);
            }
        }
        return res;
    }

    /**
     * Поиск альф для x^n. (x-a1)*(x-a2)(x-a4)(x-a8)=>x^3(a1+a2+a4+a8)
     * @param level - уровень рекурсии
     * @param degree - степень, для которой ищем альфы
     * @param curr_idx - текущий индекс в циклотомическом классе
     * @param curr_alpha - текущая альфа
     * @param curr_cycle - индекс циклотомического класса
     */
    public void findAlpha(int level, int degree, int curr_idx, int curr_alpha, int curr_cycle) {
        if (level == cycle.get(curr_cycle).size() - degree) {
            alpha.add(curr_alpha);
        } else {
            for (int i = curr_idx; i < cycle.get(curr_cycle).size(); i++) {
                findAlpha(level + 1, degree, i + 1, curr_alpha + cycle.get(curr_cycle).get(i), curr_cycle);
            }
        }
    }

    /**
     * Сумма двух полиномов по модулю 2
     * @param g1
     * @param g2
     * @return
     */
    public static boolean[] sum(boolean[] g1, boolean[] g2) {
        boolean[] res = new boolean[g1.length];
        for(int i = 0; i < g1.length; i++) {
            res[i] = (g1[i] ^ g2[i]);
        }
        return res;
    }

    /** Получение степени многочлена
     *
     * @param gx
     * @return
     */
    public static int getMaxDegree(boolean[] gx) {
        for(int i = gx.length-1; i >= 0; i--) {
            if(gx[i]) return i;
        }
        return -1;
    }

    /**
     * Перебор всех многочленов для нахождения примитивного x^m + ... + 1
     * @param i
     * @param divider
     * @param dividend
     * @return
     */
    public void BinaryBruteForce(int i, boolean[] divider, boolean[] dividend) {
        if (i == m) {
            if (dividing(dividend, divider)) {
                reducingPolies.add(divider);
            }
        } else {
            BinaryBruteForce(i + 1, divider.clone(), dividend);
            if (i != 0) divider[i] = (!divider[i]);
            BinaryBruteForce(i + 1, divider.clone(), dividend);
        }
    }
}
