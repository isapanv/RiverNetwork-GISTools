from qgis.core import (
    QgsProject,
    QgsCoordinateReferenceSystem,
    QgsPointXY,
    QgsRasterLayer,
    QgsVectorLayer,
    QgsGeometry,
    QgsField,
    QgsSpatialIndex,
    QgsFeatureRequest,
    QgsRaster,
    QgsFeature,
    QgsApplication,
    QgsPoint,
    QgsVectorLayer,
)

from qgis.analysis import QgsNativeAlgorithms
from processing_saga_nextgen.saga_nextgen_plugin import SagaNextGenAlgorithmProvider
from qgis.PyQt.QtCore import QVariant, QCoreApplication
from qgis.PyQt.QtWidgets import QAction, QFileDialog, QInputDialog, QMessageBox
import requests
import os
from pyproj import Transformer
import processing
from itertools import combinations
import random
from qgis.analysis import QgsZonalStatistics


class CustomDEMPlugin:

    def __init__(self, iface):
        self.iface = iface
        self.plugin_name = "RiverNETWORK"

    def initGui(self):
        # Для запуска плагина
        from qgis.PyQt.QtWidgets import QAction
        self.action = QAction(self.plugin_name, self.iface.mainWindow())
        self.action.triggered.connect(self.run_plugin)
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToMenu("&RiverNETWORK", self.action)

    def unload(self):
        # Удаление плагина
        self.iface.removeToolBarIcon(self.action)
        self.iface.removePluginMenu("&RiverNETWORK", self.action)

    def run_plugin(self):
        # Код плагина
        project_folder = QFileDialog.getExistingDirectory(
            None, "Выберите рабочую папку"
        )
        if not project_folder:
            QMessageBox.warning(None, "Ошибка", "Рабочая папка не выбрана. Работа плагина прекращена.")
            return

        # Создать папку "work" внутри выбранной папки
        project_folder = os.path.join(project_folder, "work")
        if not os.path.exists(project_folder):
            os.makedirs(project_folder)
        QMessageBox.information(None, "Папка установлена", f"Рабочая папка: {project_folder}")

        # Установить систему координат проекта (EPSG:3857 - Pseudo-Mercator)
        crs = QgsCoordinateReferenceSystem("EPSG:3857")
        QgsProject.instance().setCrs(crs)

        # Включить встроенные алгоритмы обработки QGIS
        QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())

        # Добавить слой OpenTopoMap
        opentopo_url = 'type=xyz&zmin=0&zmax=21&url=https://tile.opentopomap.org/{z}/{x}/{y}.png'
        opentopo_layer = QgsRasterLayer(opentopo_url, 'OpenTopoMap', 'wms')
        QgsProject.instance().addMapLayer(opentopo_layer)

        x, ok_x = QInputDialog.getDouble(None, "Координата X", "Введите координату X:", value=4316873, decimals=6)
        if not ok_x:
            QMessageBox.warning(None, "Ошибка", "Неприавильная координата X. Работа плагина прекращена.")
            return
        y, ok_y = QInputDialog.getDouble(None, "Координата Y", "Введите координату Y:", value=7711643, decimals=6)
        if not ok_y:
            QMessageBox.warning(None, "Ошибка", "Неприавильная координата Y. Работа плагина прекращена.")
            return
        x_3857, y_3857 = x, y
        QMessageBox.information(None, "Координаты установлены", f"X: {x_3857}, Y: {y_3857}")

        # Координаты озера в системе координат EPSG:3857
        # x_3857, y_3857 = 4316873, 7711643

        # Преобразовать координаты в широту и долготу (EPSG:4326)
        transformer = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)
        longitude, latitude = transformer.transform(x_3857, y_3857)

        # Определить bounding box вокруг точки с радиусом 0.5 градуса широты и долготы
        radius = 0.5
        xmin, ymin = longitude - radius, latitude - radius
        xmax, ymax = longitude + radius, latitude + radius

        # Конфигурация запроса к OpenTopography API
        api_key = 'c1fcbd0b2f691c736e3bf8c43e52a54d'
        url = (
            f"https://portal.opentopography.org/API/globaldem?"
            f"demtype=SRTMGL1"
            f"&south={ymin}&north={ymax}"
            f"&west={xmin}&east={xmax}"
            f"&outputFormat=GTiff"
            f"&API_Key={api_key}"
        )
        response = requests.get(url)

        # Сохранить скачанный DEM на диск
        output_path = f'{project_folder}srtm_output.tif'
        with open(output_path, 'wb') as f:
            f.write(response.content)

        # Добавить скачанный слой DEM в проект QGIS
        dem_layer = QgsRasterLayer(output_path, "SRTM DEM Layer")
        QgsProject.instance().addMapLayer(dem_layer)

        # Репроекция DEM из EPSG:4326 в EPSG:3857
        reprojected_relief = processing.run("gdal:warpreproject", {
            'INPUT': output_path,
            'SOURCE_CRS': QgsCoordinateReferenceSystem('EPSG:4326'),
            'TARGET_CRS': QgsCoordinateReferenceSystem('EPSG:3857'),
            'RESAMPLING': 0, 'NODATA': -9999, 'TARGET_RESOLUTION': 30,
            'OPTIONS': '', 'DATA_TYPE': 0, 'TARGET_EXTENT': None,
            'TARGET_EXTENT_CRS': None, 'MULTITHREADING': False,
            'EXTRA': '', 'OUTPUT': 'TEMPORARY_OUTPUT'})['OUTPUT']

        # Загрузка алгоритмов SAGA
        provider = SagaNextGenAlgorithmProvider()
        provider.loadAlgorithms()
        QgsApplication.processingRegistry().addProvider(provider=provider)

        # Использовать SAGA Fill Sinks для извлечения водосборов
        filled_relief = processing.run("sagang:fillsinkswangliu", {
            'ELEV': reprojected_relief, 'FILLED': 'TEMPORARY_OUTPUT',
            'FDIR': 'TEMPORARY_OUTPUT', 'WSHED': f'{project_folder}basins.sdat',
            'MINSLOPE': 0.01})['WSHED']

        # Сохранить и добавить заполненные области водосбора в проект
        basins = QgsRasterLayer(f'{project_folder}basins.sdat', 'basins')
        QgsProject.instance().addMapLayer(basins)

        # Использовать QuickOSM для запроса данных о водных путях на заданной территории
        bbox = "4261842, 4372940, 7611625, 7813231"

        # Загрузить реки
        query = processing.run('quickosm:buildqueryextent', {
            'KEY': 'waterway', 'VALUE': 'river', 'EXTENT': bbox, 'TIMEOUT': 25})
        file = processing.run("native:filedownloader", {
            'URL': query['OUTPUT_URL'], 'OUTPUT': 'TEMPORARY_OUTPUT'})['OUTPUT']
        rivers = self.iface.addVectorLayer(file + '|layername=lines', "rivers", "ogr")

        # Загрузить ручьи
        query = processing.run('quickosm:buildqueryextent', {
            'KEY': 'waterway', 'VALUE': 'stream', 'EXTENT': bbox, 'TIMEOUT': 25})
        file = processing.run("native:filedownloader", {
            'URL': query['OUTPUT_URL'], 'OUTPUT': 'TEMPORARY_OUTPUT'})['OUTPUT']
        streams = self.iface.addVectorLayer(file + '|layername=lines', "streams", "ogr")

        # Объединить слои рек и ручьев
        layer1 = QgsProject.instance().mapLayersByName("rivers")[0]
        layer2 = QgsProject.instance().mapLayersByName("streams")[0]
        merge_result = processing.run("qgis:mergevectorlayers", {
            'LAYERS': [layer1, layer2], 'CRS': layer1.crs(),
            'OUTPUT': 'TEMPORARY_OUTPUT'})['OUTPUT']
        merge_result = processing.run("native:dissolve", {'INPUT': merge_result, 'FIELD': [], 'SEPARATE_DISJOINT': False, 'OUTPUT': 'TEMPORARY_OUTPUT'})[
            'OUTPUT']
        merge_result = processing.run("native:multiparttosingleparts", {'INPUT': merge_result, 'OUTPUT': f'{project_folder}merge_result.gpkg'})['OUTPUT']

        # Добавить объединенный слой в проект
        rivers_merged = QgsVectorLayer(f'{project_folder}merge_result.gpkg', 'rivers_merged')
        QgsProject.instance().addMapLayer(rivers_merged)

        # Загрузить полигональные данные о водных объектах
        query = processing.run('quickosm:buildqueryextent', {
            'KEY': 'natural', 'VALUE': 'water', 'EXTENT': bbox, 'TIMEOUT': 25})
        file = processing.run("native:filedownloader", {
            'URL': query['OUTPUT_URL'], 'OUTPUT': f'{project_folder}water.gpkg'})['OUTPUT']
        water = self.iface.addVectorLayer(file + '|layername=multipolygons', "water", "ogr")

        # Рассчитать координаты начальных и конечных точек линий
        start_x = processing.run("native:fieldcalculator", {
            'INPUT': f'{project_folder}merge_result.gpkg', 'FIELD_NAME': 'start_x',
            'FIELD_TYPE': 0, 'FIELD_LENGTH': 0, 'FIELD_PRECISION': 0,
            'FORMULA': 'x(start_point($geometry))', 'OUTPUT': 'TEMPORARY_OUTPUT'})['OUTPUT']

        start_y = processing.run("native:fieldcalculator", {
            'INPUT': start_x, 'FIELD_NAME': 'start_y',
            'FIELD_TYPE': 0, 'FIELD_LENGTH': 0, 'FIELD_PRECISION': 0,
            'FORMULA': 'y(start_point($geometry))', 'OUTPUT': 'TEMPORARY_OUTPUT'})['OUTPUT']

        end_x = processing.run("native:fieldcalculator", {
            'INPUT': start_y, 'FIELD_NAME': 'end_x',
            'FIELD_TYPE': 0, 'FIELD_LENGTH': 0, 'FIELD_PRECISION': 0,
            'FORMULA': 'x(end_point($geometry))', 'OUTPUT': 'TEMPORARY_OUTPUT'})['OUTPUT']

        end_y = processing.run("native:fieldcalculator", {
            'INPUT': end_x, 'FIELD_NAME': 'end_y',
            'FIELD_TYPE': 0, 'FIELD_LENGTH': 0, 'FIELD_PRECISION': 0,
            'FORMULA': 'y(end_point($geometry))', 'OUTPUT': 'TEMPORARY_OUTPUT'})['OUTPUT']

        # Добавить новые поля для хранения высотных данных
        layer_provider = end_y.dataProvider()
        layer_provider.addAttributes([QgsField("start_z", QVariant.Double)])
        layer_provider.addAttributes([QgsField("end_z", QVariant.Double)])
        end_y.updateFields()

        # Начать редактирование и заполнение значений высоты
        end_y.startEditing()
        line_provider = end_y.dataProvider()

        for feature in end_y.getFeatures():
            geom = feature.geometry()
            if geom.isMultipart():
                polyline = geom.asMultiPolyline()[0]
            else:
                polyline = geom.asPolyline()

            start_point = QgsPointXY(polyline[0])
            end_point = QgsPointXY(polyline[-1])

            # Высотные данные начальной точки
            start_z = dem_layer.dataProvider().identify(start_point, QgsRaster.IdentifyFormatValue)

            start_z_idx = line_provider.fields().indexOf('start_z')
            if start_z.isValid():
                start_z_value = start_z.results()[1]
                feature['start_z'] = start_z_value

            # Высотные данные конечной точки
            end_z = dem_layer.dataProvider().identify(end_point, QgsRaster.IdentifyFormatValue)
            end_z_idx = line_provider.fields().indexOf('end_z')
            if end_z.isValid():
                end_z_value = end_z.results()[1]
                feature['end_z'] = end_z_value

            line_provider.changeAttributeValues({feature.id(): {start_z_idx: start_z_value, end_z_idx: end_z_value}})

        end_y.commitChanges()

        # Определить максимальную высоту для каждой линии
        max_z = processing.run("native:fieldcalculator", {
            'INPUT': end_y, 'FIELD_NAME': 'max_z',
            'FIELD_TYPE': 0, 'FIELD_LENGTH': 0, 'FIELD_PRECISION': 0,
            'FORMULA': 'if("start_z" > "end_z", "start_z", "end_z")',
            'OUTPUT': f'{project_folder}rivers_with_points.gpkg'})['OUTPUT']

        rivers_and_points = QgsVectorLayer(max_z, 'rivers_and_points')
        QgsProject.instance().addMapLayer(rivers_and_points)

        # Создать слой точек максимальной высоты
        point_layer = QgsVectorLayer("Point?crs=epsg:4326", "MaxHeightPoints", "memory")
        QgsProject.instance().addMapLayer(point_layer)
        layer_provider = point_layer.dataProvider()
        layer_provider.addAttributes([QgsField("x", QVariant.Double)])
        layer_provider.addAttributes([QgsField("y", QVariant.Double)])
        layer_provider.addAttributes([QgsField("z", QVariant.Double)])
        point_layer.updateFields()
        fields = point_layer.fields()

        # Получить ссылки на слои линий и точек
        line_layer_name = 'rivers_and_points'
        layers = QgsProject.instance().mapLayersByName(line_layer_name)
        layer = layers[0]

        point_layer_name = 'MaxHeightPoints'
        pointLayers = QgsProject.instance().mapLayersByName(point_layer_name)
        pointLayer = pointLayers[0]

        # Сначала собираем все конечные и начальные точки
        start_points = set()
        end_points = set()

        for feature in layer.getFeatures():
            start_x = feature['start_x']
            start_y = feature['start_y']
            end_x = feature['end_x']
            end_y = feature['end_y']

            if start_x is not None and start_y is not None:
                start_points.add((start_x, start_y))

            if end_x is not None and end_y is not None:
                end_points.add((end_x, end_y))

        pointLayer.startEditing()

        for feature in layer.getFeatures():
            if feature['max_z'] is not None:
                start_x = feature['start_x']
                start_y = feature['start_y']
                start_z = feature['start_z']
                end_x = feature['end_x']
                end_y = feature['end_y']
                end_z = feature['end_z']

                # Проверка начальной точки
                if start_x is not None and start_y is not None and start_z is not None:
                    start_point = (start_x, start_y)
                    if start_z == feature['max_z'] and start_point not in end_points:
                        point = QgsPointXY(start_x, start_y)
                        new_feature = QgsFeature()
                        new_feature.setFields(fields)
                        new_feature.setGeometry(QgsGeometry.fromPointXY(point))
                        new_feature['x'] = start_x
                        new_feature['y'] = start_y
                        new_feature['z'] = start_z
                        pointLayer.addFeature(new_feature)

                # Проверка конечной точки
                if end_x is not None and end_y is not None and end_z is not None:
                    end_point = (end_x, end_y)
                    if end_z == feature['max_z'] and end_point not in start_points:
                        point = QgsPointXY(end_x, end_y)
                        new_feature = QgsFeature()
                        new_feature.setFields(fields)
                        new_feature.setGeometry(QgsGeometry.fromPointXY(point))
                        new_feature['x'] = end_x
                        new_feature['y'] = end_y
                        new_feature['z'] = end_z
                        pointLayer.addFeature(new_feature)

        # Завершение редактирования и сохранение изменений
        pointLayer.commitChanges(True)


# Загрузка слоя точек и слоя стоимости
points_layer = QgsProject.instance().mapLayersByName('MaxHeightPoints')[0]
cost_layer = QgsProject.instance().mapLayersByName('Slope Layer')[0]

# Выбираем две случайные точки
# start_feature, end_feature = random.sample(valid_features, 2)

start_point = 'dot'
start_layer = QgsProject.instance().mapLayersByName(f'{start_point}')[0]
# end_layer = create_temp_point_layer(end_feature.geometry(), points_layer.crs())

params = {
    'INPUT_COST_RASTER': cost_layer,
    'INPUT_RASTER_BAND': 1,
    'INPUT_START_LAYER': start_layer,
    'INPUT_END_LAYER': points_layer,
    'BOOLEAN_FIND_LEAST_PATH_TO_ALL_ENDS': False,
    'BOOLEAN_OUTPUT_LINEAR_REFERENCE': False,
    'OUTPUT': 'TEMPORARY_OUTPUT'
}

# Запуск алгоритма Least-Cost Path
try:
    result = processing.run("Cost distance analysis:Least Cost Path", params)
    path_layer = result['OUTPUT']

    # Добавляем слой с результатом в проект
    QgsProject.instance().addMapLayer(path_layer)
except Exception as e:
    print(f"Ошибка при выполнении Least-Cost Path: {e}")

# удаление путей, с перепадом высот
paths_layer = QgsProject.instance().mapLayersByName('Output least cost path')[0]
elevation_layer = QgsProject.instance().mapLayersByName('SRTM DEM Layer')[0]


def calculate_minimum_elevation(raster_layer, line_geom):
    provider = raster_layer.dataProvider()
    min_elevation = float('inf')

    if line_geom.isMultipart():
        lines = line_geom.asMultiPolyline()
    else:
        lines = [line_geom.asPolyline()]

    for line in lines:
        for point in line:
            value, result = provider.sample(QgsPointXY(point.x(), point.y()), 1)  # 1 — это идентификатор первой бандовой
            if result and value is not None:
                min_elevation = min(min_elevation, value)

    return min_elevation if min_elevation != float('inf') else None


paths_to_delete = []

for path_feature in paths_layer.getFeatures():
    path_geom = path_feature.geometry()

    min_elevation = calculate_minimum_elevation(elevation_layer, path_geom)

    if min_elevation is None:
        continue

    if path_geom.isMultipart():
        first_point = path_geom.asMultiPolyline()[0][0]
        last_point = path_geom.asMultiPolyline()[-1][-1]
    else:
        first_point = path_geom.asPolyline()[0]
        last_point = path_geom.asPolyline()[-1]

    start_point = QgsPointXY(first_point.x(), first_point.y())
    end_point = QgsPointXY(last_point.x(), last_point.y())

    provider = elevation_layer.dataProvider()
    z_start, result_start = provider.sample(start_point, 1)
    z_end, result_end = provider.sample(end_point, 1)

    if not result_start or not result_end:
        continue

    z1 = min(z_start, z_end)

    if min_elevation < z1 - 15:
        paths_to_delete.append(path_feature.id())

if paths_to_delete:
    paths_layer.startEditing()
    for path_id in paths_to_delete:
        paths_layer.deleteFeature(path_id)
    paths_layer.commitChanges()

print(f"Удалено {len(paths_to_delete)} путей, не соответствующих критериям.")

# удаление путей, пересекаюжих реки
paths_layer = QgsProject.instance().mapLayersByName('Output least cost path')[0]
lines_layer = QgsProject.instance().mapLayersByName('rivers_and_points')[0]

line_index = QgsSpatialIndex(lines_layer.getFeatures())

paths_to_delete = []

for path_feature in paths_layer.getFeatures():
    path_geom = path_feature.geometry()
    path_points = path_geom.asPolyline()

    start_point = QgsGeometry.fromPointXY(path_points[0])
    end_point = QgsGeometry.fromPointXY(path_points[-1])

    candidate_ids = line_index.intersects(path_geom.boundingBox())

    for candidate_id in candidate_ids:
        line_feature = lines_layer.getFeature(candidate_id)
        line_geom = line_feature.geometry()

        if path_geom.intersects(line_geom):
            if (
                start_point.intersects(line_geom) or
                end_point.intersects(line_geom)
            ):
                continue
            paths_to_delete.append(path_feature.id())
            break

if paths_to_delete:
    paths_layer.startEditing()
    for path_id in paths_to_delete:
        paths_layer.deleteFeature(path_id)
    paths_layer.commitChanges()

print(f"Удалено {len(paths_to_delete)} пересекающихся путей (исключая общие начала/концы).")
