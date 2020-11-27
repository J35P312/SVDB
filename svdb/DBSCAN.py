import numpy

# supports up to 4D data


def x_coordinate_clustering(data, epsilon, m):
    clusters = numpy.zeros(len(data)) - 1
    cluster_id = -1
    cluster = False

    for i in range(len(data) - m + 1):
        current = data[i, :]
        points = data[i + 1:i + m, :]
        # print points
        distances = [abs(point[0] - current[0]) for point in points]
        if max(distances) < epsilon:
            if cluster:  # add to the cluster
                clusters[i + m - 1] = cluster_id
            else:  # define a new cluster
                cluster_id += 1
                cluster = True
                for j in range(i, i + m):
                    clusters[j] = cluster_id
        else:
            cluster = False
    return clusters, cluster_id


def y_coordinate_clustering(data, epsilon, m, cluster_id, clusters):
    cluster_id_list = set(clusters)
    for cluster in cluster_id_list:
        if cluster == -1:
            continue
        class_member_mask = (clusters == cluster)
        indexes = numpy.where(class_member_mask)[0]
        signals = data[class_member_mask]

        y_coordinates = [[signal[1], indexes[i]] for i, signal in enumerate(signals)]
        y_coordinates.sort(key=lambda x: x[0])

        sub_clusters = numpy.zeros(len(indexes)) - 1

        active_cluster = False
        sub_cluster_id = 0
        y_coordinates = numpy.array(y_coordinates)
        for i in range(len(y_coordinates) - m + 1):
            current = y_coordinates[i, :]
            distances = [abs(pos[0] - current[0]) for pos in y_coordinates[i + 1:i + m, :]]

            if max(distances) < epsilon:
                # add to the cluster
                if active_cluster:
                    sub_clusters[i + m - 1] = sub_cluster_id
                    # define a new cluster
                else:
                    sub_cluster_id += 1
                    active_cluster = True
                    for j in range(i, i + m):
                        sub_clusters[j] = sub_cluster_id
            else:
                active_cluster = False

        for i in range(len(sub_clusters)):
            if sub_clusters[i] == 1:
                clusters[y_coordinates[i][1]] = cluster
            elif sub_clusters[i] > -1:
                clusters[y_coordinates[i][1]] = sub_clusters[i] + cluster_id - 1
            elif sub_clusters[i] == -1:
                clusters[y_coordinates[i][1]] = -1
        if sub_cluster_id > 1:
            cluster_id += sub_cluster_id - 1
    return clusters, cluster_id


def main(data, epsilon, m):
    clusters, cluster_id = x_coordinate_clustering(data, epsilon, m)
    clusters, cluster_id = y_coordinate_clustering(data, epsilon, m, cluster_id, clusters)
    return clusters


data = numpy.array([[1, 1], [10481823, 10483880], [10481947, 10483785], [10481947, 1], [10481947, 1], [10481947, 1], [10482033, 10483984], [10482079, 10483801], [10482111, 10483972], [10482121, 10483788], [10482125, 10483769], [10482126, 10484204], [10482163, 10483811], [10482177, 10483909], [10482186, 10483906], [
                   10482191, 10483836], [10482202, 10484150], [10482262, 10483947], [10482285, 10483797], [10482342, 10483968], [10483770, 10482390], [10482390, 10483814], [10483770, 10482405], [10483769, 10482428], [10483770, 10482405], [10483770, 10482405], [10483770, 1], [10483770, 1], [10483770, 1], [10483770, 1]])
