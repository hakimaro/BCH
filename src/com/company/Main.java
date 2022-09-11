package com.company;

import java.util.Scanner;

public class Main {
    public static void main(String[] args) {
        Scanner in = new Scanner(System.in);
        System.out.print("Введите m: ");
        int m = in.nextInt();
        System.out.print("Введите d: ");
        int d = in.nextInt();
        if (d < 2 || m < 2) {
            System.out.println("Неверный ввод");
            return;
        }
        BCH test = new BCH(m, d);
        test.initReducingPolinomial();
        System.out.println("Выберите неприводимый полином");
        for(int i = 0; i < test.reducingPolies.size()/2; i++) {
            System.out.print(i + ": ");
            test.printPolynom(test.reducingPolies.get(i));
        }
        System.out.print("Выбор: ");
        test.reducingPoly = test.reducingPolies.get(in.nextInt());
        System.out.println("1. Поиск циклонических классов");
        test.initCycle();
        test.printCycle();
        System.out.println("2. Создание мультипликативной группы");
        test.initMulGroup(test.reducingPoly);
        test.printMulGroup();
        System.out.println("3. Создание порождающего многочлена");
        test.initGx();
        test.printPolynom(test.gx);
        System.out.println("--- Кодирование слова с клавиатуры ---");
        System.out.print("Введите размер слова: ");
        boolean[] s = new boolean[in.nextInt()];
        for(int i = 0; i < s.length; i++) {
            System.out.println("Введите 0 или 1");
            s[i] = in.nextInt() != 0;
        }
        boolean[] res = test.multiplying(s, test.gx);
        System.out.print("Ваше слово: ");
        test.printPolynom(s);
        System.out.print("Кодированное слово: ");
        test.printPolynom(res);
    }
}
