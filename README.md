Python-проект для анализа связности многоканальных табличных временных рядов.

Основной сценарий: взять CSV/Excel-таблицу, применить preprocessing, посчитать connectivity-метрики между каналами и сохранить HTML/Excel-отчет. Дополнительные сценарии для group comparison, HDF5/fMRI-like данных и voxel/bin-представлений остаются экспериментальными.

Создан в первую очередь для анализа временных рядов фМРТ-данных и сравнения различных метрик.

## Стабильное ядро

- загрузка табличных временных рядов из CSV и Excel;
- preprocessing: пропуски, нормализация, выбросы, AR/STL-опции;
- connectivity-метрики через единый registry;
- CLI `neweds` и Python API для одного файла;
- batch-запуск для набора табличных файлов;
- HTML/Excel-отчеты и экспорт connectivity-матриц;
- тесты для метрик, loader-ов, CLI и публичного pipeline.

Экспериментальные части вынесены ниже и не считаются основным контрактом проекта.

## Сводная таблица возможностей

### Метрики связности

| Семейство | Варианты в коде | Что возвращает | Directed | Partial/контроли | Статус |
|---|---|---|---|---|---|
| Линейная корреляция | `correlation_full` | Pearson `r` в `[-1, 1]` | нет | нет | стабильная |
| Ранговые корреляции | `correlation_spearman`, `correlation_kendall` | Spearman rho, Kendall tau-b | нет | нет | стабильные |
| Частичная корреляция | `correlation_partial` | связь после учета остальных каналов через precision matrix | нет | да | стабильная |
| Лаговая корреляция | `correlation_directed` | корреляция источника с целью на лаге | да | нет | рабочая, не стабильное ядро |
| H2 | `h2_full`, `h2_partial`, `h2_directed` | квадрат корреляционной оценки | частично | `h2_partial` | рабочие, не стабильное ядро |
| Спектральная когерентность | `coherence_full`, `coherence_partial` | magnitude-squared coherence в `[0, 1]` | нет | `coherence_partial` | `coherence_full` стабильная |
| Ordinal MI | `ordinal_full`, `ordinal_directed` | MI по порядковым паттернам Bandt-Pompe | `ordinal_directed` | нет | `ordinal_full` стабильная, directed экспериментальная |
| Mutual information | `mutinf_full`, `mutinf_partial` | kNN/KSG-подобная оценка MI | нет | `mutinf_partial` | экспериментальные |
| Distance correlation | `dcor_full`, `dcor_partial`, `dcor_directed` | dCor в `[0, 1]` | `dcor_directed` | `dcor_partial` | `dcor_full` стабильная, остальные экспериментальные |
| AH / active information | `ah_full`, `ah_partial`, `ah_directed` | направленная информационная оценка | да | `ah_partial` | экспериментальные, opt-in |
| Granger | `granger_full`, `granger_partial` | p-values F-теста Грейнджера | да | `granger_partial` | `granger_full` стабильная |
| Transfer entropy | `te_full`, `te_partial` | transfer entropy после дискретизации | да | `te_partial` | экспериментальные |

### Проверки и preprocessing

| Блок | Что делает код | Где смотреть |
|---|---|---|
| Пропуски и типы | приводит столбцы к числам, учитывает missing fraction, удаляет полностью нечисловые/NaN-ряды | `neweds.core.preprocessing`, `neweds.core.pipeline` |
| Вырожденные признаки | удаляет all-NaN, почти константные и низкодисперсные ряды | `additional_preprocessing` |
| Нормализация | z-score, robust z-score через median/MAD, ранговая нормализация | `preprocess_timeseries` |
| Выбросы | клиппинг по robust z-score, замена спайков локальной медианой | `clip_outliers_robust`, `replace_spikes_median` |
| Дрейф | линейный detrend или high-pass через вычитание скользящего среднего | `detrend_linear`, `detrend_highpass` |
| Автокорреляция и AR | оценивает `rho1`, `phi_ar1`, `tau_int`, `n_eff`, Ljung-Box; умеет AR(p)-выбеливание | `autocorr_summary`, `preprocess_timeseries` |
| Сезонность и частоты | ищет ACF/FFT-пики, может вычитать сезонную компоненту через STL | `detect_seasonality`, `fft_peaks`, `preprocess_timeseries` |
| Partial-контроли | резидуализует по контрольным переменным или использует precision matrix | `neweds.metrics.*`, `ComputationContract` |
| P-value коррекция | поддерживает `none`, `bonferroni`, `fdr_bh` для p-value-метрик | `neweds.core.pipeline` |
| Отчеты | сохраняет HTML/Excel, connectivity-матрицы и служебные сводки запуска | `neweds.reporting` |

## Установка

```bash
pip install -e .
```

Для разработки и тестов:

```bash
pip install -e ".[dev,advanced,io]"
```

`advanced` добавляет часть опциональных метрик. `io` добавляет Parquet и старый XLS через `pyarrow` и `xlrd`.

## Быстрый старт

```bash
neweds examples/demo_timeseries.csv \
  --variants correlation_full,dcor_full,ordinal_full \
  --output-dir outputs/demo
```

После запуска появится:

```text
outputs/demo/
├── report.html
├── report.xlsx
└── connectivity_exports/
```

Минимальный Python API:

```python
from neweds import AnalysisConfig, run_analysis

result = run_analysis(
    "examples/demo_timeseries.csv",
    AnalysisConfig(
        variants=["correlation_full", "dcor_full", "ordinal_full"],
    ),
)

print(result.metrics.keys())
```

## Вход и выход

Обычный вход - таблица, где строки являются временными точками, а числовые столбцы являются сигналами:

```text
signal_a, signal_b_lagged, noise_control, seasonal_component
...
```

Основной выход:

- `report.html` - интерактивный отчет;
- `report.xlsx` - табличный отчет;
- экспортированные connectivity-матрицы;
- служебные сводки запуска.

## Экспериментальные сценарии

- `neweds-group` для case/control-сравнения;
- HDF5/fMRI-like эксперименты;
- выравнивание voxel/bin-пространства;
- validation-сценарии.

Эти сценарии полезны для раннего анализа и проверки pipeline на исследовательских данных, но их результаты нужно валидировать отдельно. Подробнее: [Ограничения](docs/limitations.md).

## Документация

- [Архитектура](docs/architecture.md)
- [Метрики](docs/metrics.md)
- [Математика и предварительный анализ данных](docs/math_and_data_analysis.md)
- [Group pipeline](docs/group_pipeline.md)
- [Ограничения](docs/limitations.md)
- [Refactoring story](docs/refactoring_story.md)
- [Demo](examples/README.md)
