import torch
import torch.nn as nn
import torch.nn.functional as F

class CustomModel(nn.Module):
    def __init__(self):
        super(CustomModel, self).__init__()
        self.conv1 = nn.Conv2d(in_channels=1, out_channels=10, kernel_size=(5, 5), padding='same')
        self.bn1 = nn.BatchNorm2d(10)
        self.conv2 = nn.Conv2d(10, 20, kernel_size=(7, 7), padding='same')
        self.bn2 = nn.BatchNorm2d(20)
        self.pool = nn.MaxPool2d(kernel_size=(3, 3), stride=(1, 1), padding=0)
        self.dropout = nn.Dropout(0.5)
        self.lstm = nn.LSTM(input_size=900, hidden_size=8, batch_first=True, bidirectional=True)
        self.fc1 = nn.Linear(16, 10)  # 16 units total (8 in each direction)
        self.fc2 = nn.Linear(10, 1)
        self.fc3 = nn.Linear(1, 1)

    def forward(self, x):
        #print('initial: ', x.shape)
        x = x.view(-1, 1, 30, 47)  # Reshape input
        x = F.relu(self.bn1(self.conv1(x)))
        #print('conv1: ', x.shape)
        x = F.relu(self.bn2(self.conv2(x)))
        #print('conv2: ', x.shape)
        x = self.pool(x)
        #print('maxpooling: ', x.shape)
        x = x.permute(0, 2, 1, 3)  # Flatten for LSTM
        x = x.reshape(-1, 28, 900)
        #print('flatten:', x.shape)

        x = self.dropout(x)
        x, _ = self.lstm(x)  # Reshape for LSTM
        #print('lstm:', x.shape)
        x = self.dropout(x[:, -1, :])  # Use output of last LSTM sequence
        x = F.relu(self.fc1(x))
        x = F.softplus(self.fc2(x))
        x = torch.sigmoid(self.fc3(x))
        #print('output:', x.shape)
        return x