'''
ANN model
'''
import torch
import torch.nn as nn

class PositionWiseFeedForward(nn.Module):
    
    def __init__(self, num_hidden, num_ff):
        super(PositionWiseFeedForward, self).__init__()
        
        self.W_in = nn.Linear(num_hidden, num_ff, bias=True)
        self.W_out = nn.Linear(num_ff, num_hidden, bias=True)
        #activation function: gaussian error linear unit
        self.act = torch.nn.GELU()

    def forward(self, h_V):
        h = self.act(self.W_in(h_V))
        h = self.W_out(h)
        return h
