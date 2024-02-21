import pandas as pd
import numpy as np
import cv2 as cv


from collections import namedtuple
from time import perf_counter

error = 25
backup_iterations = 0

## POINTS IN THE PHANTOM
df_org = pd.read_csv('data/reference.csv', index_col=0)
original_x, original_y = df_org.iloc[:, 0], df_org.iloc[:, 1]
original_points = np.array([original_x, original_y]).T


image_org = cv.imread('data/sample.tif')
image_org = cv.cvtColor(image_org, cv.COLOR_BGR2GRAY)

X, Y = np.meshgrid(np.arange(-512, 512), np.arange(-512, 512))
image_org[np.sqrt(X**2 + Y**2) > 979 / 2] = 0

image = cv.medianBlur(image_org, 7)

kernel = np.array([
    [1, 0, -1],
    [3, 0, -3],
    [1, 0, -1]
])
image = cv.GaussianBlur(image, (7, 7), 0)

image1 = cv.filter2D(image, -1, kernel)
image2 = cv.filter2D(image, -1, -kernel)
image3 = cv.filter2D(image, -1, kernel.T)
image4 = cv.filter2D(image, -1, -(kernel.T))

image = image1 + image2 + image3 + image4
 
image_org[image > 50] = 0

circles = cv.HoughCircles(image_org, cv.HOUGH_GRADIENT, 1, 20, None, 100, 4, 6, 8)

detected_points = []
show_image = cv.cvtColor(image, cv.COLOR_GRAY2RGB)
if circles is not None:
    circles = np.uint16(np.around(circles))
    for i in circles[0, :]:
        detected_points.append((i[0], i[1]))
        cv.circle(show_image, (i[0], i[1]), i[2], (0, 0, 255), 1)
else:
    print("No circles detected")
    exit()

detected_points = np.array(detected_points, dtype='float32')
detected_x, detected_y = detected_points[:, 0], detected_points[:, 1]

with open('data/reference.txt', 'w') as f:
    f.write(f'{len(original_points)}\n')
    np.savetxt(f, original_points)

with open('data/sample.txt', 'w') as f:
    f.write(f'{len(detected_points)}\n')
    np.savetxt(f, detected_points)


cv.imshow('Processed image', cv.resize(show_image, (512, 512)))
cv.waitKey(0)
cv.destroyAllWindows()

Match = namedtuple('Match', ['original', 'detected', 'size'])

def expandMatch(match, new, step, original_distance_matrix, detected_distance_matrix):
    for i, j in zip(match.detected, match.original):
        if np.linalg.norm(detected_distance_matrix[i, step] - original_distance_matrix[j, new]) > error:
            return False
    return True

def getDistanceMatrix(points):
    return np.array([np.linalg.norm(p - points, axis=1) for p in points])


print(f'''
    number of detected points: {len(detected_points)}
    number of original points: {len(original_points)}
    ratio of detected points: {len(detected_points) / len(original_points)}
''')

original_distance_matrix = getDistanceMatrix(original_points)
detected_distance_matrix = getDistanceMatrix(detected_points)

n, m = len(detected_points), len(original_points)

matches = []
for i in range(m):
    for j in range(m):
        if abs(detected_distance_matrix[0, 1] - original_distance_matrix[i, j]) < error:
            matches.append(Match(original=(i, j), detected=(0, 1), size=2))

candidatesInEachStep = [len(matches)]
biggestCandidateInEachStep = [max(match.size for match in matches)]

print('Number of matches in each step:')
print('\t 1. ', len(matches), biggestCandidateInEachStep[0])

start_time = perf_counter()
for step in range(2, len(detected_points)):
    ###### Update Section ######
    new_matches = []
    for match in matches:
        foundMatch = False
        for i in range(m):
            if i not in match and expandMatch(match, i, step, original_distance_matrix, detected_distance_matrix):
                new_matches.append(Match(
                    original=(*match.original, i),
                    detected=(*match.detected, step),
                    size=match.size + 1
                ))
                foundMatch = True
        if not foundMatch:
            new_matches.append(match)

    matches = new_matches

    ###### Filter Section ######
    biggestMatch = max(matches, key=lambda match: match.size)

    matches = [
            match for match in matches if match.detected[-1] > step - backup_iterations
        ] + [
            match for match in matches if match.size == biggestMatch.size and match.detected[-1] <= step - backup_iterations
        ]
    ###########################


    ## count matches of each size  
    candidatesInEachStep.append(len(matches))
    biggestCandidateInEachStep.append(biggestMatch.size)


end_time = perf_counter()
print()
print(f'Number of matches: {len(matches)}')
if len(matches) != 0:
    print(f'Biggest candidate: {biggestCandidateInEachStep[-1]}')


bestMatches = [match for match in matches if match.size == biggestCandidateInEachStep[-1]]
print(f'Number of best matches: {len(bestMatches)}')
print(f'Time: {end_time - start_time}')


m = len(detected_points)
ground_truth = np.array(bestMatches[0].original) * m + np.array(bestMatches[0].detected)
ground_truth.sort()

with open('data/ground_truth.txt', 'w') as f:
    f.write(f'{len(ground_truth)}\n')
    np.savetxt(f, ground_truth)

