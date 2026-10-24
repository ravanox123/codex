# DNA Lab — проект по биоинформатике

Небольшой Python-проект для анализа последовательностей ДНК в формате FASTA. Работает на Python 3.10+ без внешних зависимостей.

Возможности:

- Чтение нескольких последовательностей FASTA и проверка символов IUPAC.
- Длина, количество A/C/G/T, число неоднозначных оснований и GC-состав.
- Обратная комплементарная последовательность.
- Перевод ДНК в белок со стандартным генетическим кодом и выбором рамки.
- Поиск точного мотива, включая перекрывающиеся совпадения.
- Сохранение отчёта в JSON.

## Запуск

Из папки проекта:

```bash
python3 -m dna_lab.cli examples/sequences.fasta --motif ATG
python3 -m dna_lab.cli examples/sequences.fasta --frame 1 --output report.json
```

Дополнительно можно установить CLI: `python3 -m pip install -e .`, затем использовать `dna-lab` вместо `python3 -m dna_lab.cli`.

## Проверка

```bash
python3 -m unittest discover -s tests -v
```

## Как трактовать результаты

GC-процент рассчитывается среди A/C/G/T; неоднозначные основания исключаются из знаменателя. Если A/C/G/T отсутствуют, результат `null`. Позиции мотивов начинаются с 1; символы IUPAC в мотиве сопоставляются буквально, без расширения вариантов. Рамка перевода — смещение 0, 1 или 2 от начала последовательности. Стоп-кодон обозначается `*`, неоднозначный кодон — `X`; неполный последний кодон пропускается. Перевод продолжается после стоп-кодонов и не выполняет поиск генов или ORF.

## Структура

```text
dna_lab/core.py       # алгоритмы анализа
dna_lab/cli.py        # командная строка
examples/            # демонстрационные FASTA
tests/               # проверки алгоритмов и чтения файлов
```

## WES Lab desktop app

A local Python application for WES variant review, SFARI gene overlap, graphs,
coverage and report exports is now included. Start with
`./scripts/start_wes_lab.sh --example`. See [WES_LAB.md](WES_LAB.md) for supported
inputs, the sequential BAM/annotation workflow and interpretation limits.
